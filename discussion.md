### Introduction

**SMILES** (simplified molecular-input line-entry system) is a classic notation used in organic chemistry to describe molecular formulas and uniform scientific communication. FSMILES notation is developped by Feng at al. [https://arxiv.org/pdf/2305.10133.pdf] in the context of Ligvo3DModel algorithm development allowing generation of new 3D molecules.

Example of SMILES and FSMILES notations for **2-{4-\[(4'-Methoxybiphenyl-3-Yl)sulfonyl]piperazin-1-Yl}-3-(4-Methoxyphenyl)pyrazine**

**2D-structure:**

<img src = "2-{4-[(4'-Methoxybiphenyl-3-Yl)sulfonyl]piperazin-1-Yl}-3-(4-Methoxyphenyl)pyrazine.png" />

*(image source: https://pubchem.ncbi.nlm.nih.gov/compound/56955953#section=Structures)*

**SMILES:**
*COC1=CC=C(C=C1)C2=CC(=CC=C2)S(=O)(=O)N3CCN(CC3)C4=NC=CN=C4C5=CC=C(C=C5)OC*

**FSMILES:**
*'start_0'O_0=_0S_0(\[\*\])_0(\[\*\])_0=_0O_0'sep_0’N_61_0C_6C_6N_6(\[\*\])_0C_6C_61_0'sep_0’c_61_0n_6c_6c_6n_6c_61_0\[\*\]_0'sep_0’c_61_0c_6c_6c_6([\*])_0c_6c_61_0'sep_0’O_0C_0'sep_0’c_61_0c_6c_6c_6c_6(\[\*\])_0c_61_0'sep_0’c_61_0c_6c_6c_6(\[\*\])_0c_6c_61_0'sep_ 0’O_0C_0
'sep_0''end_0'*

-------

### FSMILES features and advantages

**FSMILES** (fragment-based simplified molecular-input line-entry system) notation has several key features that distinguish it from SMILES:

* All ring tokens as well as acyclic fragments can be clearly distinguished due to separator **'sep_0'**;

* Each chemical element contains information about **ring size** (number of elements in corresponding ring tokens). Acyclic fragments have ring size equal to zero;

* The **connection points** are explicitely denoted by ([\*]) (or [\*] at if a connection point is associated with the last element of the chemical group);

* The **aromaticity** of each token is explicitely defined (elements in uppercase for saturated or acyclic tokens, elements in lowercase for aromatic tokens);

* FSMILES notation starts from acyclic fragment located in the middle of molecule (if exists).

This makes FSMILES notation longer and less readable for human but useful for **Lingvo3DModel** and potentially for ML and DL algorithms.

Some advatages highlighted by the authors include:

* FSMILES integrates **local** spherical **coordinate system.** As bond lengths and bond angles are rigid, and bond angles can be computed knowing atom type and ring size, FSMILES notation can be converted in **global** Euclidean *coordinate system,** so it is indirectly integrated in FSMILES notation. This combination of local and global coordinate systems helps predict more accurately chemical substructures of generated 3D molecules.

* When generating new molecules with Lingvo3DModel, FSMILES helps prioritize **ring closure** with accurate bond angles and reasonable ring size.

--------

### FSMILES for transformers

FSMILES has some potential advantages when used with transformers. Transformers are typically used with text information. The first step of input data transformation is **tokenization** when the text is divided into tokens (in classic transformers each token is a word). Then, tokens are processed by two parallel transformations: sense embedding and positional embedding (keeps information about position of token in the data).

Chemical formulas in FSMILES notation are more suitable to use with transformers then formulas in SMILES notaion: chemical groups can be considered as tokens, separators between chemical groups allow token distinction. After tokenization, sense embedding vector gets information about local coordinates and positional embedding vector -- about global coordinates of a token respectively.
