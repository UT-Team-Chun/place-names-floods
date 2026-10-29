# Toponymic audit of the readings: coding guide (amendment B-12, item B21 and addendum 1)

Purpose: an independent check, by specialists in the languages and the place-name literature concerned, of how the study
read the names of settlements. This guide replaces `v68_audit_guide.md` (amendment B-11). The sample now also contains
names that the study did not read as water-related, so that the audit can estimate how many water-related names the study
missed, and it adds names from the countries that carry most of the study's estimates near settlements.

Please do not consult the study's code repository, the paper, the lexicon or the study's readings until your coding is
complete. Please also do not consult reports, maps or images of floods at these places. The study's readings of this
sample are kept in a separate file that is not public and will be compared with your codes only after the audit.

## The sample

`v69_audit_sample.csv` lists 300 settlements in a random order: `item` (1 to 300), `gid` (GeoNames identifier), `name`
(the GeoNames name, with diacritics), `asciiname` (the ASCII form the study read), `cc` (country) and the coordinates.
Some of the names were read by the study as containing a water-related element and some were not; the file does not say
which. Most names come from Cambodia, Myanmar, Pakistan, Bangladesh, Australia, Germany, Great Britain and Portugal, with
others from about twenty further countries.

## What to record for each row

Code every element of the name, not only the name as a whole.

1. `language` — the language (or language layer) of the name and of each element (for example Khmer, Burmese, Urdu,
   Punjabi, Saraiki, Sindhi, Bengali, Low German, Irish, English as transferred in Australia).
2. `elements` — the elements of the name, each with its meaning and language, from an attested early form where one
   exists (give the form and its date if you use one).
3. `water_related` — for each element, whether it refers to water or wet ground (river, stream, lake, spring, well, marsh,
   fen, bog, reed bed, river island, water meadow, floodplain, bank, dike, polder, ford, bridge, ferry, landing, valley,
   hollow, low ground): `yes`, `no` or `uncertain`.
4. `class` — for water-related elements: `water`, `wet` (wetland), `riverine` (river island, water meadow, floodplain,
   bank, dike, polder), `crossing`, `low` (valley, hollow, lowland).
5. `qualifier` — whether an element is a relative qualifier that distinguishes one of a pair or group of settlements
   (such as English *lower*, *nether*, *upper*, German *nieder*, *ober*, Slovene *spodnji*, Finnish *ala-*): `yes` or `no`.
   Code a qualifier as water-related only if you judge that it refers to low or wet ground rather than to the position of
   one settlement relative to another.
6. `site_or_feature` — for each water-related element: `site` when the element describes the ground on which the original
   settlement stood (for example "the settlement on the river island"); `feature` when it names a distinct feature beside
   or near the settlement (a ford, a brook, a marsh next to a settlement on a dry rise); `both` when the name supports
   either reading; `unknown` when the sources do not say. Use the sources rather than the modern map.
7. `motivation` — the kind of name: `topographical` (describes the landscape), `habitative` (a word for a farm, village or
   town), `eponymous` (from a person), `commemorative` or `transferred` (from another place), or `other`.
8. `homograph` — whether the string could also be another word or a personal name, title or transferred name (for
   example German *-au* from Slavic *-ov*, Persian *-ab* in personal names): `yes` with a note, or `no`.
9. `source` — the dictionary or survey consulted (author, title, page or entry). Please do not rely on Wikipedia alone.
10. `confidence` — `high`, `medium` or `low`.
11. `comment` — optional.

## Coders, double coding and disagreements

- At least two coders per language block, and at least 20 % of the names in each language coded independently by two
  coders. The pairs of codes are kept as they are, and Cohen's kappa is reported by language.
- Where the two codes disagree, the coders discuss the name and record an agreed code with a short reason; where they
  cannot agree, a third specialist decides. Both original codes are kept.
- Coders work from this guide and the sources; they do not see each other's codes until both have finished.

## Standard sources (examples)

- England: Watts, *The Cambridge Dictionary of English Place-Names* (2004); Smith, *English Place-Name Elements* (1956);
  the English Place-Name Society county surveys; the Key to English Place-Names.
- German: Niemeyer, *Deutsches Ortsnamenbuch* (2012); the Bavarian and Austrian historical place-name books.
- France: Nègre, *Toponymie générale de la France*; Dauzat and Rostaing, *Dictionnaire étymologique des noms de lieux en
  France*.
- Slovenia: Snoj, *Etimološki slovar slovenskih zemljepisnih imen* (2009). Poland: Rymut (ed.), *Nazwy miejscowe Polski*.
- Sweden: *Svenskt ortnamnslexikon*. Finland: *Suomalainen paikannimikirja*.
- Cambodia, Myanmar, Bangladesh and Pakistan: specialists in Khmer, Burmese, Bengali, Punjabi, Saraiki and Sindhi place
  names; national gazetteers and district gazetteers.

## How the codes will be used

Your codes will be compared with the study's readings to estimate, by form, language and country: the share of names the
study read as water-related whose matched element you code as water-related (precision), the share of names the study did
not read as water-related that contain a water-related element (the rate of missed names), and how often the matched
element is a qualifier, a homograph, or names a feature rather than the settlement's site. The study's corrections for
misreading will be recomputed with these rates, or labelled illustrative where the audit does not cover a form.
