# Parse-Nextflow-Integrated-Pipeline-CRIANN

Adaptation de **ParseNIP** (Parse Nextflow-Integrated Pipeline) pour tourner sur le cluster CRIANN (partition `austral`), sans Docker.

Pipeline original : [StevenWingett/ParseNIP](https://github.com/StevenWingett/ParseNIP)
Construit split-pipe (Parse Biosciences) pour l'analyse scRNA-seq (Parse Evercode) : construction d'index génomique, mapping, QC, objets Scanpy/Seurat, démultiplexage par sample.
---

## Contenu du dépôt

| Fichier | Rôle |
|---|---|
| `main.nf` | Pipeline Nextflow principal (ParseNIP) |
| `nextflow_austral.config` | Configuration adaptée à CRIANN (sans Docker/Apptainer, environnements conda) |
| `ParseNIP.sh` | Script de soumission SLURM (v1.5) |
| `bin/check_setup.py` | Vérifie la cohérence des paramètres/fichiers avant de lancer split-pipe (index/mapping uniquement — voir limite connue ci-dessous) |
| `bin/fastq_sep_groups_v1_5.py` | Script Parse Biosciences : sépare les reads par sample well group (démultiplexage) |
| `bin/concat_fastq.sh` | Fusionne les FASTQ démultiplexés par sample à travers les sublibrairies |
| `modules/local/splitpipe/index/main.nf` | Process Nextflow : construction d'index split-pipe via `--fasta`/`--genes` |
| `modules/local/splitpipe/index_gfasta/main.nf` | Process Nextflow : construction d'index split-pipe via `--gfasta` (sans GTF) |
| `modules/local/python3_check_setup/check_setup_build_index/main.nf` | Vérification de setup pour le mode `--fasta`/`--gtf` |
| `modules/local/python3_check_setup/check_setup_build_index_gfasta/main.nf` | Vérification de setup pour le mode `--gfasta` |
| `modules/local/parse_demux/demultiplex/main.nf` | Process Nextflow : démultiplexage par sample d'une sublibrairie (`fastq_sep_groups_v1_5.py`) |
| `modules/local/parse_demux/concatenate_demux/main.nf` | Process Nextflow : fusion des FASTQ démultiplexés par sample, toutes sublibrairies confondues |
| `Edit_fna_gtf_TOGA2.sh` | Prépare FASTA+GTF issus de TOGA2 pour split-pipe, mode `--fasta`/`--gtf` (voir ci-dessous) |
| `Edit_fna_gfasta_TOGA2.sh` | Prépare un FASTA seul (transcrits épissés TOGA2) pour split-pipe, mode `--gfasta` (voir ci-dessous) |
| `fix_gfasta_headers.sh` | Corrige a posteriori les headers d'un FASTA `--gfasta` déjà généré (voir Dépannage) |
| `get_mitochondrion.py` | Récupère et formate le génome mitochondrial depuis NCBI |

---

## Pré-requis sur CRIANN

Pas de Docker disponible sur le cluster : le pipeline tourne avec des environnements **conda** activés à la volée dans chaque tâche SLURM, plus quelques modules CRIANN pour l'orchestrateur.

### Modules
```bash
module load py_env/miniconda3/25.7.0 tools/nextflow/25.10.4
```

### Environnements conda requis

| Environnement | Contenu | Utilisé pour |
|---|---|---|
| `spipe_env` | split-pipe 1.4.0, fastqc, trim-galore, scanpy 1.11.5, seaborn, biopython, numpy, pandas | FastQC, Trim Galore, construction d'index, mapping, combine, QC split-pipe, Scanpy, démultiplexage |
| `seurat_env` | Seurat 5.0.0 (R) | Construction des objets Seurat |
| `multiqc_env` | MultiQC | Rapport final agrégé |

L'activation de chaque environnement est gérée automatiquement par `nextflow_austral.config` (directive `beforeScript`) — aucune activation manuelle n'est nécessaire dans le script de soumission SLURM.

---

## Deux modes de construction d'index : `--fasta`/`--gtf` vs `--gfasta`

Le pipeline supporte deux façons de construire l'index génomique split-pipe, **mutuellement exclusives** :

- **Mode `--fasta`/`--gtf`** (historique) : génome + annotation classique, coordonnées GTF exprimées sur le référentiel du FASTA fourni.
- **Mode `--gfasta`** (ajouté pour les génomes TOGA2) : un seul FASTA, sans GTF. split-pipe génère automatiquement un enregistrement gene/exon par séquence FASTA, couvrant sa longueur totale. **Nécessaire quand le FASTA contient des séquences déjà épissées/extraites par transcrit** (comme `nucleotide.fa` produit par TOGA2), car dans ce cas un GTF en coordonnées génomiques (scaffold complet, avec introns) ne correspond plus à la longueur des séquences et casse silencieusement la construction de l'index STAR (exons "hors chromosome", silencieusement ignorés).

Utiliser le mode `--fasta`/`--gtf` avec un FASTA de séquences pré-extraites est le piège classique qui a motivé l'ajout du mode `--gfasta` à ce fork.

**⚠️ Limite connue du mode `--gfasta` avec TOGA2** : `nucleotide.fa` de TOGA2 ne contient que le **CDS** (séquence codante), pas l'ARNm complet (UTR exclues — TOGA ne transfère que les annotations protéine-codantes). Comme la chimie Parse Evercode capture les transcrits par leur queue poly-A, la majorité des reads séquencent la région **3' UTR**. Une référence sans UTR peut donc donner un taux de mapping très bas (observé : ~8%). Vérifier si le run TOGA2 a été fait avec l'option d'annotation UTR (nouvelle en TOGA2) avant de construire l'index, sous peine de perdre l'essentiel du signal en mapping.

---

## Démultiplexage par sample (`--demultiplex`)

Étape indépendante du flux index/mapping habituel : sépare les FASTQ **bruts** (jamais trimmés ni concaténés) par sample well group, via le script Parse Biosciences `fastq_sep_groups_v1_5.py`.

- Tourne une fois par sublibrairie → `<outdir>/demultiplexing/<sublibrary_id>/`
- Puis fusionne les FASTQ par sample à travers toutes les sublibrairies (adapté de `concat_fastq.sh`) → `<outdir>/demultiplexing/combined/`
- Réutilise directement `--samp_list` comme fichier de groupes (`-g/--gfile` de `fastq_sep_groups_v1_5.py`) : le format `sample_id wells` est le même que celui déjà validé par `check_setup.py`.
- Kit par défaut : `WT` (paramètre `--kit`, voir `fastq_sep_groups_v1_5.py --kit_list` pour les valeurs valides).
- Peut tourner seul (juste `--fastq`/`--samp_list`/`--demultiplex`, sans construction d'index ni mapping) ou en complément d'un run `index`/`mapping`/`all` classique — les deux flux ne dépendent pas l'un de l'autre : fournir `--fastq`/`--samp_list` sans index (ni existant ni à construire) ne déclenche **plus** le mapping (voir Dépannage), à condition que `--demultiplex` soit actif.

```bash
# Démultiplexage seul
nextflow run main.nf -config nextflow_austral.config -profile austral \
  --fastq fastqs/ --samp_list samples.tsv --demultiplex --kit WT --outdir results/

# En complément d'un run complet
nextflow run main.nf -config nextflow_austral.config -profile austral \
  --fasta genome.fa --gtf genome.gtf --genome_name MonGenome \
  --fastq fastqs/ --samp_list samples.tsv --demultiplex \
  --outdir results/
```

**⚠️ Limite connue** : en run "démultiplexage seul" (sans index), `check_setup.py` n'est **pas** appelé pour valider `--fastq`/`--samp_list` (le pipeline part directement sur `CHECK_SETUP_DUMMY_RUN`, faute d'index à vérifier). La seule validation en amont vient de `fastq_sep_groups_v1_5.py` lui-même. Si besoin d'une vraie vérification FASTQ/samplesheet pour ce cas, il faudrait étendre `check_setup.py` avec un 4ème mode dédié (non fait à ce jour).

---

## Utilisation de `main.nf`

### 1. Construire un index génomique split-pipe uniquement

**Mode `--fasta`/`--gtf` :**
```bash
nextflow run main.nf \
  -config nextflow_austral.config \
  -profile austral \
  --fasta genome.fa --gtf genome.gtf --genome_name MonGenome \
  --outdir results/ \
  -resume
```

**Mode `--gfasta` :**
```bash
nextflow run main.nf \
  -config nextflow_austral.config \
  -profile austral \
  --gfasta genome_forGfasta.fa --genome_name MonGenome \
  --outdir results/ \
  -resume
```

### 2. Mapper sur un index déjà construit
```bash
nextflow run main.nf \
  -config nextflow_austral.config \
  -profile austral \
  --genome_dir results/index_MonGenome/ \
  --fastq fastqs/ --samp_list samples.tsv \
  --outdir results/ \
  -resume
```
(identique quel que soit le mode utilisé pour construire l'index — `SPLITPIPE_MAP`/`SPLITPIPE_COMBINE` sont agnostiques du mode.)

### 3. Construire l'index et mapper dans le même run
```bash
nextflow run main.nf \
  -config nextflow_austral.config \
  -profile austral \
  --fasta genome.fa --gtf genome.gtf --genome_name MonGenome \
  --fastq fastqs/ --samp_list samples.tsv \
  --trim \
  --outdir results/ \
  -resume
```
(remplacer `--fasta genome.fa --gtf genome.gtf` par `--gfasta genome_forGfasta.fa` pour le mode gfasta)

Voir `nextflow run main.nf --help` pour la liste complète des paramètres.

---

## Soumission via SLURM

```bash
sbatch ParseNIP.sh -r RUN -g REFGENOME -s STEP [-m gtf|gfasta] [-o outputdir] [-d] [-k kit]
```

- `-r, --run` : nom du lot de FASTQ à traiter
- `-g, --genome` : nom du génome de référence à utiliser (**obligatoire pour `-s index/mapping/all`, pas pour `-s demultiplex`**)
- `-s, --step` : `index` (construction d'index seule), `mapping` (mapping seul sur index existant), `all` (pipeline complet : index, mapping, Seurat), ou `demultiplex` (démultiplexage seul, sans toucher au génome — implique `-d` automatiquement)
- `-m, --mode` : `gtf` (défaut, `--fasta`/`--gtf`) ou `gfasta` (`--gfasta`, sans GTF)
- `-o, --outdir` : nom du dossier de sortie (défaut : `results`)
- `-d, --demultiplex` : active le démultiplexage **en plus** de l'étape choisie via `-s` (combinable avec `index`/`mapping`/`all`)
- `-k, --kit` : nom du kit pour le démultiplexage (défaut : `WT`)

Consulter l'aide complète (avec exemples) : `bash ParseNIP.sh -h` — **toujours via `bash`, jamais `sbatch`** : `sbatch` soumet le job à la file d'attente au lieu d'exécuter le script immédiatement, donc le message d'aide finirait dans les fichiers `.out`/`.err` au lieu du terminal, après un délai d'attente inutile.

**Exemples :**
```bash
# Tout, avec démultiplexage (index + mapping + Seurat + démultiplexage par sample)
sbatch ParseNIP.sh -r Test -g Myotis_myotis -s all -m gfasta -d

# Démultiplexage seul, pas de génome nécessaire
sbatch ParseNIP.sh -r Test -s demultiplex -k WT_mega

# Tout, sans démultiplexage (comportement classique)
sbatch ParseNIP.sh -r Test -g Myotis_myotis -s all -m gfasta

# Mapping seul, index déjà construit
sbatch ParseNIP.sh -r Test -g Myotis_myotis -s mapping
```

Le script construit automatiquement les chemins attendus selon le mode, et **valide strictement qu'exactement un fichier correspond à chaque motif** avant de lancer Nextflow (`resolve_one_file`, depuis v1.3) — échec explicite et immédiat si aucun ou plusieurs fichiers matchent, plutôt qu'un échec tardif et confus côté Nextflow :
- mode `gtf` : `${REFGENOMEPATH}/${REFGENOME}_*.fa.gz` + `${REFGENOMEPATH}/${REFGENOME}_*.gtf.gz`
- mode `gfasta` : `${REFGENOMEPATH}/TOGA2/${REFGENOME}*_forGfasta.fa.gz`

Cette résolution ne s'exécute que pour `-s index`/`all` (pas pour `mapping`/`demultiplex`, qui n'en ont pas besoin).

### Dossiers de sortie

Deux dossiers de destination distincts, selon ce qui a été produit :

| Run | `${FOLDER}/OutputSplitpipe/${RUN}/` | `${FOLDER}/OutputDemultiplexing/${RUN}/` |
|---|---|---|
| `-s demultiplex` | rien | `results/demultiplexing/...` |
| `-s index`/`mapping`/`all` (sans `-d`) | `results/...` (tout) | rien |
| `-s all -d` (ou `index`/`mapping` `-d`) | `results/...` **sans** le sous-dossier `demultiplexing/` | `demultiplexing/...` (extrait, pas dupliqué) |

Quand `-d` est combiné à une autre étape, le sous-dossier `demultiplexing/` est **déplacé** vers `OutputDemultiplexing`, pas copié aux deux endroits.

---

## Préparation des génomes de référence (TOGA2 → split-pipe)

TOGA2 fournit deux sorties utilisables comme point de départ :
- `nucleotide.fa.gz` : séquences **transcrits épissés (CDS uniquement)** — courtes, sans intron, sans UTR
- `query_annotation.gtf.gz` : annotation en **coordonnées génomiques/scaffold** (avec introns) — beaucoup plus longues

Ces deux fichiers ne sont **pas directement compatibles** entre eux (systèmes de coordonnées différents) : c'est pourquoi il existe deux scripts de préparation selon la source utilisée.

### `Edit_fna_gtf_TOGA2.sh` (mode `--fasta`/`--gtf`, génomes complets en coordonnées scaffold)

```bash
bash Edit_fna_gtf_TOGA2.sh -g GENOME -mt ACCESSION
```

- `-g, --genome` : préfixe du génome (ex. `Myotis_myotis_Bat1k`), attend en entrée `GENOME.gtf` et `GENOME_nucleotide.fa` **en coordonnées scaffold cohérentes entre les deux fichiers**
- `-mt, --mtaccession` : numéro d'accession NCBI du génome mitochondrial complet de l'espèce (ex. `KT901455`)

Étapes : nettoyage des sorties précédentes, harmonisation `#`→`_` dans les attributs GTF et les en-têtes FASTA, report du `transcript_id` sur les lignes gene/exon/CDS, synthèse de lignes `gene` si absentes, ajout du génome mitochondrial (`get_mitochondrion.py`), vérification de cohérence FASTA/GTF (9 champs, noms de séquences), compression.

**Fichiers produits** : `GENOME_nucleotide_filtered_mtDNA_forParse.fa.gz` + `GENOME_filtered_mtDNA_forParse.gtf.gz` → à passer à `main.nf` via `--fasta`/`--gtf`.

### `Edit_fna_gfasta_TOGA2.sh` (mode `--gfasta`, transcrits épissés TOGA2)

```bash
bash Edit_fna_gfasta_TOGA2.sh -g GENOME -mt ACCESSION
```

- `-g, --genome` : préfixe du génome (ex. `Myotis_myotis_Bat1k`), attend en entrée `GENOME_nucleotide.fa` (le `nucleotide.fa.gz` de TOGA2, transcrits épissés)
- `-mt, --mtaccession` : identique à ci-dessus

Étapes : **nettoyage des en-têtes par liste blanche** — tout caractère hors `[A-Za-z0-9._-]` (virgule, `$`, parenthèses, guillemets, espaces...) est remplacé par `_`, un par un, sans supprimer ni fusionner de champ. C'est important : les headers TOGA2 ont la forme `ENScode#GeneName#nombre[,nombre...]`, où le champ numérique final **distingue de vrais paralogues** (jusqu'à plusieurs loci par gène) — le supprimer plutôt que le nettoyer recrée des doublons de header (observé : ~1700 collisions sur un seul génome). Suivent : suppression des séquences vides, **détection et arrêt sur en-têtes dupliqués** (collision de nom de chromosome dans l'index STAR sinon), ajout du génome mitochondrial (mêmes gènes que ci-dessus, sans le GTF), compression. Aucun GTF n'est manipulé ni produit — `--gfasta` en génère un automatiquement en interne à partir de la longueur de chaque séquence.

**Fichier produit** : `GENOME_nucleotide_forGfasta.fa.gz` → à passer à `main.nf` via `--gfasta`.

### `get_mitochondrion.py`

```bash
python3 get_mitochondrion.py ACCESSION SPECIES_NAME
```

- Télécharge (ou réutilise si déjà présent en local) l'enregistrement GenBank de l'accession donnée
- Extrait les features `CDS`, `tRNA`, `rRNA`
- Désambiguïse automatiquement les gènes dupliqués (isoaccepteurs tRNA-Leu/Ser) en `Leu1`/`Leu2`, `Ser1`/`Ser2`, etc.
- Produit un FASTA (`SPECIES_ACCESSION_genes.fa`) et un GTF (`SPECIES_ACCESSION_genes.gtf`) formatés correctement pour split-pipe
  (le GTF n'est utilisé que par `Edit_fna_gtf_TOGA2.sh`, ignoré par `Edit_fna_gfasta_TOGA2.sh`)

**Dépendance** : nécessite `biopython`, présent dans `spipe_env`.

---

## Dépannage

**`Dataframe load failed; Error tokenizing data`** lors de `split-pipe --mode mkref` : le GTF contient des lignes avec plus de 9 colonnes (tabulations dans le champ attributs, souvent après une édition manuelle). Vérifier avec :
```bash
awk -F'\t' '{count[NF]++} END {for (n in count) printf "Number of lines with %d fields: %d\n", n, count[n]}' fichier.gtf
```

**Avertissements STAR `exon end = X is larger than the chromosome ... length = Y, will skip this exon`** (des milliers de lignes) lors de `--runMode genomeGenerate` : signe que le FASTA et le GTF fournis via `--fasta`/`--gtf` ne partagent pas le même système de coordonnées — typiquement un FASTA de séquences déjà extraites/épissées (courtes) combiné à un GTF en coordonnées génomiques complètes (grandes). C'est exactement le cas des sorties brutes TOGA2 (`nucleotide.fa.gz` + `query_annotation.gtf.gz`) : **ne pas les combiner directement**. Utiliser `Edit_fna_gfasta_TOGA2.sh` + le mode `--gfasta` à la place.

**`ValueError: invalid literal for int() with base 10: '2300"'`** (ou similaire) dans `molinfo.py` / `join_csv_chunk_files`, en toute fin de `SPLITPIPE_MAP` (après un run presque complet) : un `gene_name` dérivé d'un header `--gfasta` contient un caractère qui casse le parsing CSV naïf interne de split-pipe (typiquement une virgule héritée du champ numérique TOGA2 `ENScode#GeneName#nombre,nombre`). Corriger les headers du FASTA déjà produit avec `fix_gfasta_headers_v3.sh` (liste blanche `[A-Za-z0-9._-]`, ne supprime aucun champ distinctif — voir avertissement ci-dessus sur les paralogues), puis **reconstruire l'index** avant de relancer le mapping. Pour repartir de zéro proprement, ce nettoyage est aussi intégré dans `Edit_fna_gfasta_TOGA2.sh`.

**Taux de mapping très bas (à un chiffre %) en mode `--gfasta`** avec un génome TOGA2 : voir l'avertissement dans la section "Deux modes de construction d'index" ci-dessus — `nucleotide.fa` de TOGA2 est probablement CDS-only, sans UTR, alors que Parse Evercode échantillonne majoritairement la région 3' UTR.

**`ERROR ~ Argument of 'file()' function cannot be null` (ligne du bloc `CHECK_SETUP_PREBUILT_INDEX` dans `main.nf`)** lors d'un run `--demultiplex` seul (sans `--genome_dir` ni construction d'index) : bug corrigé — `checkParameters()` déclenchait `perform_mapping = true` dès que `--fastq`/`--samp_list` étaient fournis, même sans index disponible, puis tentait `file(params.genome_dir)` avec `genome_dir` resté `null`. Depuis correctif, `perform_mapping` ne se déclenche que si `--genome_dir` **ou** une construction d'index (`--fasta`/`--gtf`/`--gfasta`) est aussi présent ; sinon, avec `--demultiplex` actif, le run se limite bien au démultiplexage seul.

**`CONCATENATE_DEMULTIPLEXED` termine sans erreur mais `combined/` reste vide** : dû à un `ln -s` pointant vers un glob qui n'a rien matché — `ln -s` ne vérifie jamais que sa cible existe, donc l'échec est silencieux (pas de crash, juste un lien symbolique cassé). Corrigé : le process utilise maintenant `find` récursif (au lieu d'un glob à un seul niveau) et échoue explicitement (`exit 1` + message clair) si aucun fichier `.fastq.gz` n'est trouvé, ou si `concat_fastq.sh` tourne sans rien produire.

---

## Auteurs

Nicolas Nesi — Université Caen Normandie, DYNAMICURE INSERM UMR 1311
