# Etymology audit: coding guide (amendment B-6, item 5)

Purpose: an independent check of the place-name readings used in the study. You will read the Wikipedia
article(s) listed for each settlement and record what the article says about the origin of the name.
Please do not consult the study's code repository, the paper or the study's own codes
(`v44_etym_codes.csv`) until your coding is complete; the codes to be compared are kept in a separate file.

For each row of `v51_audit_sample.csv`:

1. Open the article(s) in the `articles` column and find the section on the name (etymology, toponymy, name).
2. `status`:
   - `coded` — the article gives a derivation of the name (even if tentative);
   - `no-derivation` — the article has no statement on where the name comes from;
   - `unclear` — the article discusses the name but the derivation cannot be classified.
3. `categories` (only for `coded`; several allowed, separated by `;`). Code what the recorded derivation
   refers to, not what the modern place is like:
   - `water` — a river, stream, brook, lake, pond, spring or well;
   - `wet` — a marsh, fen, bog, swamp, mud, reeds, wet meadow;
   - `riverine` — an island in a river, a river meadow, floodplain, riverbank, sandbank, dyke or polder;
   - `crossing` — a ford, bridge, ferry, landing place or river port;
   - `low` — a valley, lowland, hollow, basin, or "lower";
   - `high` — a hill, mountain, height, rock, or "upper";
   - `other` — anything else (a person, saint, tribe, building, plant other than reeds, animal, etc.).
   A name taken from the name of a river counts as `water`.
4. `quote`: copy the shortest sentence(s) from the article that support the codes.
5. `comment`: optional (for example, competing derivations or a derivation you consider doubtful).

The sample is 100 names drawn at random within strata of the study's design (seed 20261028); the strata and
the study's codes are in `v51_audit_key.csv` and are used only after your coding is complete. Agreement will
be reported as Cohen's kappa for `status` and for each category.
