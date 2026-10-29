# Toponymic audit of the readings: coding guide (amendment B-11, item B11)

Purpose: an independent check, by specialists in the languages and the place-name literature concerned, of how the study
read the names of settlements. This guide replaces the guide of amendment B-6 (`v51_audit_guide.md`), which asked coders
to record what Wikipedia says. Here we ask you to judge each name against standard sources.

Please do not consult the study's code repository, the paper, the lexicon or the study's readings until your coding is
complete. The study's readings of this sample are kept in a separate file that is not public and will be compared with
your codes only after the audit.

## The sample

`v68_audit_sample.csv` lists 150 settlements: `gid` (GeoNames identifier), `name` (the GeoNames name, with diacritics),
`asciiname` (the ASCII form the study read), `cc` (country), and the coordinates. The names were drawn within strata of
matched form and country, weighted by how much each form contributes to the study's results (the share of cells nearest a
name read as water-related and the share of flooded cells within 500 m of such names), up to three names per stratum.
Most names come from Australia, Germany, Great Britain, Pakistan, Portugal and Bangladesh, with a few from twelve other
countries.

## What to record for each row

1. `language` — the language (or language layer) of the name or of each element, if known (for example Urdu, Punjabi,
   Saraiki, Sindhi, Bengali, Low German, Bavarian, Irish, English as transferred in Australia).
2. `elements` — the elements of the name, each with its meaning and language, from an attested early form where one
   exists (give the form and its date if you use one).
3. `water_related` — for each element, whether it refers to water or wet ground (river, stream, lake, spring, well,
   marsh, fen, bog, reed bed, river island, water meadow, floodplain, bank, dike, polder, ford, bridge, ferry, landing,
   valley, hollow, low ground): `yes`, `no` or `uncertain`.
4. `class` — for water-related elements: `water`, `wet` (wetland), `riverine` (river island, water meadow, floodplain,
   bank, dike, polder), `crossing`, `low` (valley, hollow, lowland, or a relative qualifier such as "lower").
5. `site_or_feature` — whether the element names the settlement's own site or a feature beside it or nearby
   (for example a marsh next to a settlement on a dry rise): `site`, `feature`, `both` or `unknown`.
6. `homograph` — whether the matched string could also be another word or a personal name, title or transferred name
   (for example German *-au* from Slavic *-ov*, Persian *-ab* in personal names, English names transferred to Australia):
   `yes` with a note, or `no`.
7. `source` — the dictionary or survey consulted (author, title, page or entry). Please do not rely on Wikipedia alone.
8. `confidence` — `high`, `medium` or `low`.
9. `comment` — optional.

## Standard sources (examples)

- England: Watts, *The Cambridge Dictionary of English Place-Names* (2004); the English Place-Name Society county surveys;
  the Key to English Place-Names.
- German: Niemeyer, *Deutsches Ortsnamenbuch* (2012); the Bavarian and Austrian historical place-name books.
- France: Nègre, *Toponymie générale de la France*; Dauzat and Rostaing, *Dictionnaire étymologique des noms de lieux en
  France*.
- Slovenia: Snoj, *Etimološki slovar slovenskih zemljepisnih imen* (2009).
- Poland: Rymut (ed.), *Nazwy miejscowe Polski*.
- Sweden: *Svenskt ortnamnslexikon*. Finland: *Suomalainen paikannimikirja*.
- Bangladesh and Pakistan: specialists in Bengali, Punjabi, Saraiki and Sindhi place names; district gazetteers.

## How the codes will be used

Your codes will be compared with the study's readings to estimate, by form and by country, how often the element that the
study matched is water-related in your judgement (precision), with Cohen's kappa between coders where two coders code the
same names. This sample cannot estimate how many water-related names the study missed (recall). The study's corrections
for misreading will be recomputed with these rates, or labelled illustrative where the audit does not cover a form.
