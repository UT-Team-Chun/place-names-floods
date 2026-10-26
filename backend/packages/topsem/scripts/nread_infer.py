"""N-read の推論（Miyabi-G・vLLM）。入力の名前表 → 語義カテゴリと確信度。

入力 parquet: 列 gid, asciiname, name, cc, country（一意な P 名）。
出力: ``<out>/nread.parquet``（gid, category, confidence, parsed_ok, raw）と ``nread_meta.json``
（モデル snapshot・プロンプト sha256・版・件数・処理時間）。

    python nread_infer.py --names data/nread/names.parquet --out runs/<job> \
        --model google/gemma-4-31B-it --revision <snapshot>
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from topsem.world.nread import MODEL, REVISION, SYSTEM, parse, prompt_sha256, user_prompt


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--names", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", default=MODEL)
    ap.add_argument("--revision", default=REVISION, help="HF snapshot のコミット（凍結値）")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--max-model-len", type=int, default=1024)
    args = ap.parse_args()

    from vllm import LLM, SamplingParams

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    names = pd.read_parquet(args.names)
    if args.limit:
        names = names.head(args.limit)
    convs = [
        [{"role": "system", "content": SYSTEM},
         {"role": "user", "content": user_prompt(str(r.asciiname), str(r.name), str(r.country))}]
        for r in names.itertuples(index=False)
    ]
    t0 = time.time()
    llm = LLM(model=args.model, revision=args.revision, tokenizer_revision=args.revision,
              dtype="bfloat16", max_model_len=args.max_model_len,
              gpu_memory_utilization=0.90, enable_prefix_caching=True,
              limit_mm_per_prompt={"image": 0, "audio": 0})
    t_load = time.time() - t0
    sp = SamplingParams(temperature=0.0, max_tokens=48, seed=0)
    res = llm.chat(convs, sp, use_tqdm=False)
    texts = [r.outputs[0].text for r in res]
    parsed = [parse(t) for t in texts]
    df = pd.DataFrame({
        "gid": names["gid"].to_numpy(),
        "category": [p[0] for p in parsed],
        "confidence": [p[1] for p in parsed],
        "parsed_ok": [p[2] for p in parsed],
        "raw": texts,
    })
    df.to_parquet(out / "nread.parquet", index=False)
    meta = {
        "model": args.model,
        "revision": args.revision,
        "prompt_sha256": prompt_sha256(),
        "n_names": len(df),
        "parsed_ok_rate": float(df["parsed_ok"].mean()) if len(df) else None,
        "category_counts": df["category"].value_counts().to_dict(),
        "load_sec": round(t_load, 1),
        "infer_sec": round(time.time() - t0 - t_load, 1),
        "sampling": {"temperature": 0.0, "max_tokens": 48, "seed": 0},
    }
    try:
        import transformers
        import vllm
        meta["versions"] = {"vllm": vllm.__version__, "transformers": transformers.__version__}
    except ImportError:
        pass
    (out / "nread_meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
