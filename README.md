# Protein-Modeling-of-ponA-geneac
Protein-Modeling-of-ponA-geneac

#Project Details
The objective was to identify a non-homologous essential gene of unique metabolism pathways, which can be used as potential drug target against Haemophilus Influenzae. The genome of Haemophilus influenzae Rd KW20 retrieved from NCBI. The genome was analyzed using different bioinformatics tools & techniques.

The shortlisted genes of Haemophilus influenzae Rd KW20 were shortlisted in substractive genome analysis methods.
The protein structure of short listed genes were designed using Modeler 9.1. To get the 3D protein structure, the programming language python was implanted.


# implementation of python in the project
MODELLER is used for homology or comparative modeling of protein three-dimensional structures of protein molecules. The user provides an alignment of a sequence to be modeled with known related structures and MODELLER automatically calculates a model containing all non-hydrogen atoms. MODELLER implements comparative protein structure modeling by satisfaction of spatial restraints (3,4), and can perform many additional tasks, including de novo modeling of loops in protein structures, optimization of various models of protein structure with respect to a flexibly defined objective function, multiple alignment of protein sequences and/or structures, clustering, searching of sequence databases, comparison of protein structures, etc. MODELLER is available for most Unix/Linux systems, Windows, and Mac.


The best tamplet of pon-A gene is selected and the most predicated structure of the gene is developed by implementing Python programming languge.


#Where and How Python is used for the project
Basically there are five type of modelling like:

# Basic Modeling

1. Searching for structures related to TvLDH

First, it is necessary to put the target TvLDH sequence into the PIR format readable by MODELLER (file "TvLDH.ali").

>P1;TvLDH
sequence:TvLDH:::::::0.00: 0.00
MSEAAHVLITGAAGQIGYILSHWIASGELYGDRQVYLHLLDIPPAMNRLTALTMELEDCAFPHLAGFVATTDPKA
AFKDIDCAFLVASMPLKPGQVRADLISSNSVIFKNTGEYLSKWAKPSVKVLVIGNPDNTNCEIAMLHAKNLKPEN
FSSLSMLDQNRAYYEVASKLGVDVKDVHDIIVWGNHGESMVADLTQATFTKEGKTQKVVDVLDHDYVFDTFFKKI
GHRAWDILEHRGFTSAASPTKAAIQHMKAWLFGTAPGEVLSMGIPVPEGNPYGIKPGVVFSFPCNVDKEGKIHVV
EGFKVNDWLREKLDFTEKDLFHEKEIALNHLAQGG*
File: TvLDH.ali

The first line contains the sequence code, in the format ">P1;code". The second line with ten fields separated by colons generally contains information about the structure file, if applicable. Only two of these fields are used for sequences, "sequence" (indicating that the file contains a sequence without known structure) and "TvLDH" (the model file name). The rest of the file contains the sequence of TvLDH, with "*" marking its end. The standard one-letter amino acid codes are used. (Note that they must be upper case; some lower case letters are used for non-standard residues. See the file modlib/restyp.lib in the Modeller distribution for more information.)

A search for potentially related sequences of known structure can be performed by the profile.build() command of MODELLER. The following script, taken line by line, does the following (see file "build_profile.py"):

Initializes the 'environment' for this modeling run, by creating a new 'environ' object. Almost all MODELLER scripts require this step, as the new object (which we call here 'env', but you can call it anything you like) is needed to build most other useful objects.
Creates a new 'sequence_db' object, calling it 'sdb'. 'sequence_db' objects are used to contain large databases of protein sequences.
Reads a text format file containing non-redundant PDB sequences at 95% sequence identity into the sdb database. The sequences can be found in the file "pdb_95.pir" (which can be downloaded using the link at the top of this page). Like the previously-created alignment, this file is in PIR format. Sequences which have fewer than 30 or more than 4000 residues are discarded, and non-standard residues are removed.
Writes a binary machine-specific file containing all sequences read in the previous step.
Reads the binary format file back in. Note that if you plan to use the same database several times, you should use the previous two steps only the first time, to produce the binary database. On subsequent runs, you can omit those two steps and use the binary file directly, since reading the binary file is a lot faster than reading the PIR file.
Creates a new 'alignment' object, calling it 'aln', reads our query sequence "TvLDH" from the file "TvLDH.ali", and converts it to a profile 'prf'. Profiles contain similar information to alignments, but are more compact and better for sequence database searching.
Searches the sequence database 'sdb' for our query profile 'prf'. Matches from the sequence database are added to the profile.
Writes a profile of the query sequence and its homologs (see file "build_profile.prf"). The equivalent information is also written out in standard alignment format.
from modeller import *

log.verbose()
env = environ()

#-- Prepare the input files

#-- Read in the sequence database
sdb = sequence_db(env)
sdb.read(seq_database_file='pdb_95.pir', seq_database_format='PIR',
         chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

#-- Write the sequence database in binary form
sdb.write(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
          chains_list='ALL')

#-- Now, read in the binary database
sdb.read(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
         chains_list='ALL')

#-- Read in the target sequence/alignment
aln = alignment(env)
aln.append(file='TvLDH.ali', alignment_format='PIR', align_codes='ALL')

#-- Convert the input sequence/alignment into
#   profile format
prf = aln.to_profile()

#-- Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

#-- Write out the profile in text format
prf.write(file='build_profile.prf', profile_format='TEXT')

#-- Convert the profile back to alignment format
aln = prf.to_alignment()

#-- Write out the alignment file
aln.write(file='build_profile.ali', alignment_format='PIR')

File: build_profile.py

Note that while this script is written in the Python programming language, it uses Modeller-specific commands, and should therefore be run by using a command similar to the following at your command line:

mod9v1 build_profile.py

Note that the "mod9v1" script runs Modeller version 9v1. If you using a different version of Modeller, you will need to adjust this command accordingly - for example, if you have Modeller version 9.10 installed, use "mod9.10" instead.

(You can get a command line using xterm or GNOME Terminal in Linux, Terminal in Mac OS X, or the 'Modeller' link from your Start Menu in Windows. For more information on running Modeller, see the release notes. For more information on using Python, see the Python web site. Note that you can use other Python modules within your Modeller scripts, if Python is correctly installed on your system.)

The profile.build() command has many options. In this example rr_file is set to use the BLOSUM62 similarity matrix (file "blosum62.sim.mat" provided in the MODELLER distribution). Accordingly, the parameters matrix_offset and gap_penalties_1d are set to the appropriate values for the BLOSUM62 matrix. For this example, we will run only one search iteration by setting the parameter n_prof_iterations equal to 1. Thus, there is no need for checking the profile for deviation (check_profile set to False). Finally, the parameter max_aln_evalue is set to 0.01, indicating that only sequences with e-values smaller than or equal to 0.01 will be included in the final profile.

2. Selecting a template

The output of the "build_profile.py" script is written to the "build_profile.log" file. MODELLER always produces a log file. Errors and warnings in log files can be found by searching for the "_E>" and "_W>" strings, respectively. MODELLER also writes the profile in text format to the "build_profile.prf" file. An extract (omitting the aligned sequences) of the output file can be seen next. The first 6 commented lines indicate the input parameters used in MODELLER to build the profile. Subsequent lines correspond to the detected similarities by profile.build().

# Number of sequences:     30
# Length of profile  :    335
# N_PROF_ITERATIONS  :      1
# GAP_PENALTIES_1D   :   -900.0   -50.0
# MATRIX_OFFSET      :    0.0
# RR_FILE            : ${MODINSTALL8v1}/modlib//as1.sim.mat
    1 TvLDH                                    S     0   335     1   335     0     0     0    0.    0.0
    2 1a5z                                     X     1   312    75   242    63   229   164   28.   0.83E-08
    3 1b8pA                                    X     1   327     7   331     6   325   316   42.    0.0
    4 1bdmA                                    X     1   318     1   325     1   310   309   45.    0.0
    5 1t2dA                                    X     1   315     5   256     4   250   238   25.   0.66E-04
    6 1civA                                    X     1   374     6   334    33   358   325   35.    0.0
    7 2cmd                                     X     1   312     7   320     3   303   289   27.   0.16E-05
    8 1o6zA                                    X     1   303     7   320     3   287   278   26.   0.27E-05
    9 1ur5A                                    X     1   299    13   191     9   171   158   31.   0.25E-02
   10 1guzA                                    X     1   305    13   301     8   280   265   25.   0.28E-08
   11 1gv0A                                    X     1   301    13   323     8   289   274   26.   0.28E-04
   12 1hyeA                                    X     1   307     7   191     3   183   173   29.   0.14E-07
   13 1i0zA                                    X     1   332    85   300    94   304   207   25.   0.66E-05
   14 1i10A                                    X     1   331    85   295    93   298   196   26.   0.86E-05
   15 1ldnA                                    X     1   316    78   298    73   301   214   26.   0.19E-03
   16 6ldh                                     X     1   329    47   301    56   302   244   23.   0.17E-02
   17 2ldx                                     X     1   331    66   306    67   306   227   26.   0.25E-04
   18 5ldh                                     X     1   333    85   300    94   304   207   26.   0.30E-05
   19 9ldtA                                    X     1   331    85   301    93   304   207   26.   0.10E-05
   20 1llc                                     X     1   321    64   239    53   234   164   26.   0.20E-03
   21 1lldA                                    X     1   313    13   242     9   233   216   31.   0.31E-07
   22 5mdhA                                    X     1   333     2   332     1   331   328   44.    0.0
   23 7mdhA                                    X     1   351     6   334    14   339   325   34.    0.0
   24 1mldA                                    X     1   313     5   198     1   189   183   26.   0.13E-05
   25 1oc4A                                    X     1   315     5   191     4   186   174   28.   0.18E-04
   26 1ojuA                                    X     1   294    78   320    68   285   218   28.   0.43E-05
   27 1pzgA                                    X     1   327    74   191    71   190   114   30.   0.16E-06
   28 1smkA                                    X     1   313     7   202     4   198   188   34.    0.0
   29 1sovA                                    X     1   316    81   256    76   248   160   27.   0.93E-03
   30 1y6jA                                    X     1   289    77   191    58   167   109   33.   0.32E-05 
File: build_profile.prf

The most important columns in the profile.build() output are the second, tenth, eleventh and twelfth columns. The second column reports the code of the PDB sequence that was compared with the target sequence. The PDB code in each line is the representative of a group of PDB sequences that share 95% or more sequence identity to each other and have less than 30 residues or 30% sequence length difference. The eleventh column reports the percentage sequence identities between TvLDH and a PDB sequence normalized by the lengths of the alignment (indicated in the tenth column). In general, a sequence identity value above approximately 25% indicates a potential template unless the alignment is short (i.e., less than 100 residues). A better measure of the significance of the alignment is given in the twelfth column by the e-value of the alignment. In this example, six PDB sequences show very significant similarities to the query sequence with e-values equal to 0. As expected, all the hits correspond to malate dehydrogenases (1bdm:A, 5mdh:A, 1b8p:A, 1civ:A, 7mdh:A, and 1smk:A). To select the most appropriate template for our query sequence over the six similar structures, we will use the alignment.compare_structures() command to assess the structural and sequence similarity between the possible templates (file "compare.py").

from modeller import *

env = environ()
aln = alignment(env)
for (pdb, chain) in (('1b8p', 'A'), ('1bdm', 'A'), ('1civ', 'A'),
                     ('5mdh', 'A'), ('7mdh', 'A'), ('1smk', 'A')):
    m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
aln.malign()
aln.malign3d()
aln.compare_structures()
aln.id_table(matrix_file='family.mat')
env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)
File: compare.py

In this case, we create an (initially empty) alignment object 'aln' and then use a Python 'for' loop to instruct MODELLER to read each of the PDB files. (Note that in order for this to work, you must have all of the PDB files in the same directory as this script, either downloaded from the PDB website or from the archive linked at the top of this page.) We use the model_segment argument to ask only for a single chain to be read from each PDB file (by default, all chains are read from the file). As each structure is read in, we use the append_model method to add the structure to the alignment.

At the end of the loop, all of the structures are in the alignment, but they are not ideally aligned to each other (append_model creates a simple 1:1 alignment with no gaps). Therefore, we improve this alignment by using malign to calculate a multiple sequence alignment. The malign3d command then performs an iterative least-squares superposition of the six 3D structures, using the multiple sequence alignment as its starting point. The compare_structures command compares the structures according to the alignment constructed by malign3d. It does not make an alignment, but it calculates the RMS and DRMS deviations between atomic positions and distances, differences between the mainchain and sidechain dihedral angles, percentage sequence identities, and several other measures. Finally, the id_table command writes a file with pairwise sequence distances that can be used directly as the input to the dendrogram command (or the clustering programs in the PHYLIP package). dendrogram calculates a clustering tree from the input matrix of pairwise distances, which helps visualizing differences among the template candidates. Excerpts from the log file are shown below (file "compare.log").

Sequence identity comparison (ID_TABLE):

   Diagonal       ... number of residues;
   Upper triangle ... number of identical residues;
   Lower triangle ... % sequence identity, id/min(length).

         1b8pA @11bdmA @11civA @25mdhA @27mdhA @21smkA @2
1b8pA @1      327     194     147     151     153      49
1bdmA @1       61     318     152     167     155      56
1civA @2       45      48     374     139     304      53
5mdhA @2       46      53      42     333     139      57
7mdhA @2       47      49      87      42     351      48
1smkA @2       16      18      17      18      15     313


Weighted pair-group average clustering based on a distance matrix:


                                           .----------------------- 1b8pA @1.9    39.0000
                                           |
                                  .-------------------------------- 1bdmA @1.8    50.5000
                                  |
                              .------------------------------------ 5mdhA @2.4    55.3750
                              |
                              |                                .--- 1civA @2.8    13.0000
                              |                                |
        .---------------------------------------------------------- 7mdhA @2.4    83.2500
        |
      .------------------------------------------------------------ 1smkA @2.5

      +----+----+----+----+----+----+----+----+----+----+----+----+
    86.0600   73.4150   60.7700   48.1250   35.4800   22.8350   10.1900
         79.7375   67.0925   54.4475   41.8025   29.1575   16.5125
                                                                                                            
Excerpts of the file compare.log

The comparison above shows that 1civ:A and 7mdh:A are almost identical, both sequentially and structurally. However, 7mdh:A has a better crystallographic resolution (2.4Å versus 2.8Å), eliminating 1civ:A. A second group of structures (5mdh:A, 1bdm:A, and 1b8p:A) share some similarities. From this group, 5mdh:A has the poorest resolution leaving for consideration only 1bdm:A and 1b8p:A. 1smk:A is the most diverse structure of the whole set of possible templates. However, it is the one with the lowest sequence identity (34%) to the query sequence. We finally pick 1bdm:A over 1b8p:A and 7mdh:A because of its better crystallographic R-factor (16.9%) and higher overall sequence identity to the query sequence (45%).

3. Aligning TvLDF with the template

A good way of aligning the sequence of TvLDH with the structure of 1bdm:A is the align2d() command in MODELLER. Although align2d() is based on a dynamic programming algorithm, it is different from standard sequence-sequence alignment methods because it takes into account structural information from the template when constructing an alignment. This task is achieved through a variable gap penalty function that tends to place gaps in solvent exposed and curved regions, outside secondary structure segments, and between two positions that are close in space. As a result, the alignment errors are reduced by approximately one third relative to those that occur with standard sequence alignment techniques. This improvement becomes more important as the similarity between the sequences decreases and the number of gaps increases. In the current example, the template-target similarity is so high that almost any alignment method with reasonable parameters will result in the same alignment. The following MODELLER script aligns the TvLDH sequence in file "TvLDH.ali" with the 1bdm:A structure in the PDB file "1bdm.pdb" (file "align2d.py").

from modeller import *

env = environ()
aln = alignment(env)
mdl = model(env, file='1bdm', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='1bdmA', atom_files='1bdm.pdb')
aln.append(file='TvLDH.ali', align_codes='TvLDH')
aln.align2d()
aln.write(file='TvLDH-1bdmA.ali', alignment_format='PIR')
aln.write(file='TvLDH-1bdmA.pap', alignment_format='PAP')
File: align2d.py

In this script, we again create an 'environ' object to use as input to later commands. We create an empty alignment 'aln', and then a new protein model 'mdl', into which we read the chain A segment of the 1bdm PDB structure file. The append_model() command transfers the PDB sequence of this model to the alignment and assigns it the name of "1bdmA" (align_codes). Then we add the "TvLDH" sequence from file "TvLDH.seq" to the alignment, using the append() command. The align2d() command is then executed to align the two sequences. Finally, the alignment is written out in two formats, PIR ("TvLDH-1bdmA.ali") and PAP ("TvLDH-1bdmA.pap"). The PIR format is used by MODELLER in the subsequent model building stage, while the PAP alignment format is easier to inspect visually. Due to the high target-template similarity, there are only a few gaps in the alignment. In the PAP format, all identical positions are marked with a "*" (file "TvLDH-1bdmA.pap").

 _aln.pos         10        20        30        40        50        60
1bdmA     MKAPVRVAVTGAAGQIGYSLLFRIAAGEMLGKDQPVILQLLEIPQAMKALEGVVMELEDCAFPLLAGL 
TvLDH     MSEAAHVLITGAAGQIGYILSHWIASGELYG-DRQVYLHLLDIPPAMNRLTALTMELEDCAFPHLAGF 
 _consrvd *     *  ********* *   ** **  * *  * * ** ** **  *    ********* ***

 _aln.p   70        80        90       100       110       120       130
1bdmA     EATDDPDVAFKDADYALLVGAAPRL---------QVNGKIFTEQGRALAEVAKKDVKVLVVGNPANTN 
TvLDH     VATTDPKAAFKDIDCAFLVASMPLKPGQVRADLISSNSVIFKNTGEYLSKWAKPSVKVLVIGNPDNTN 
 _consrvd  ** **  **** * * **   *             *  **   *  *   **  ***** *** ***

 _aln.pos  140       150       160       170       180       190       200
1bdmA     ALIAYKNAPGLNPRNFTAMTRLDHNRAKAQLAKKTGTGVDRIRRMTVWGNHSSIMFPDLFHAEVD--- 
TvLDH     CEIAMLHAKNLKPENFSSLSMLDQNRAYYEVASKLGVDVKDVHDIIVWGNHGESMVADLTQATFTKEG 
 _consrvd   **   *  * * **     ** ***    * * *  *       *****   *  **  *

 _aln.pos    210       220       230       240       250       260       270
1bdmA     -GRPALELVDMEWYEKVFIPTVAQRGAAIIQARGASSAASAANAAIEHIRDWALGTPEGDWVSMAVPS 
TvLDH     KTQKVVDVLDHDYVFDTFFKKIGHRAWDILEHRGFTSAASPTKAAIQHMKAWLFGTAPGEVLSMGIPV 
 _consrvd          *       *      *   *   **  ****   *** *   *  **  *   **  *

 _aln.pos      280       290       300       310       320       330
1bdmA     Q--GEYGIPEGIVYSFPVTAK-DGAYRVVEGLEINEFARKRMEITAQELLDEMEQVKAL--GLI 
TvLDH     PEGNPYGIKPGVVFSFPCNVDKEGKIHVVEGFKVNDWLREKLDFTEKDLFHEKEIALNHLAQGG 
 _consrvd      ***  * * ***      *   ****   *   *     *   *  * *
File: TvLDH-1bdmA.pap

4. Model building

Once a target-template alignment is constructed, MODELLER calculates a 3D model of the target completely automatically, using its automodel class. The following script will generate five similar models of TvLDH based on the 1bdm:A template structure and the alignment in file "TvLDH-1bdmA.ali" (file "model-single.py").

from modeller import *
from modeller.automodel import *
#from modeller import soap_protein_od

env = environ()
a = automodel(env, alnfile='TvLDH-1bdmA.ali',
              knowns='1bdmA', sequence='TvLDH',
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
File: model-single.py

The first line loads in the automodel class and prepares it for use. We then create an automodel object, call it 'a', and set parameters to guide the model building procedure. alnfile names the file that contains the target-template alignment in the PIR format. knowns defines the known template structure(s) in alnfile ("TvLDH-1bdmA.ali"). sequence defines the name of the target sequence in alnfile. assess_methods requests one or more assessment scores (discussed in more detail in the next section). starting_model and ending_model define the number of models that are calculated (their indices will run from 1 to 5). The last line in the file calls the make method that actually calculates the models.

The most important output file is "model-single.log", which reports warnings, errors and other useful information including the input restraints used for modeling that remain violated in the final model. The last few lines from this log file are shown below.

>> Summary of successfully produced models:
Filename                          molpdf     DOPE score    GA341 score
----------------------------------------------------------------------
TvLDH.B99990001.pdb           1763.56104   -38079.76172        1.00000
TvLDH.B99990002.pdb           1560.93396   -38515.98047        1.00000
TvLDH.B99990003.pdb           1712.44104   -37984.30859        1.00000
TvLDH.B99990004.pdb           1720.70801   -37869.91406        1.00000
TvLDH.B99990005.pdb           1840.91772   -38052.00781        1.00000
Excerpts of the file model-single.log

As you can see, the log file gives a summary of all the models built. For each model, it lists the file name, which contains the coordinates of the model in PDB format. The models can be viewed by any program that reads the PDB format, such as Chimera. The log also shows the score(s) of each model, which are further discussed below. (Note that the actual numbers may be slightly different on your machine - this is nothing to worry about.)

5. Model evaluation

If several models are calculated for the same target, the "best" model can be selected in several ways. For example, you could pick the model with the lowest value of the MODELLER objective function or the DOPE or SOAP assessment scores, or with the highest GA341 assessment score, which are reported at the end of the log file, above. (The objective function, molpdf, is always calculated, and is also reported in a REMARK in each generated PDB file. The DOPE, SOAP, and GA341 scores, or any other assessment scores, are only calculated if you list them in assess_methods. To calculate the SOAP score, you will first need to download the SOAP-Protein potential file from the SOAP website, then uncomment the SOAP-related lines in model-single.py by removing the '#' characters.) The molpdf, DOPE, and SOAP scores are not 'absolute' measures, in the sense that they can only be used to rank models calculated from the same alignment. Other scores are transferable. For example GA341 scores always range from 0.0 (worst) to 1.0 (native-like); however GA341 is not as good as DOPE or SOAP at distinguishing 'good' models from 'bad' models.

Once a final model is selected, it can be further assessed in many ways. Links to programs for model assessment can be found in the MODEL EVALUATION section on this page.

Before any external evaluation of the model, one should check the log file from the modeling run for runtime errors ("model-single.log") and restraint violations (see the MODELLER manual for details). The file "evaluate_model.py" evaluates an input model with the DOPE potential. (Note that here we arbitrarily picked the second generated model - you may want to try other models.)

from modeller import *
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# read model file
mdl = complete_pdb(env, 'TvLDH.B99990002.pdb')

# Assess with DOPE:
s = selection(mdl)   # all atom selection
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='TvLDH.profile',
              normalize_profile=True, smoothing_window=15)
File: evaluate_model.py

In this script we use the complete_pdb script to read in a PDB file and prepare it for energy calculations (this automatically allows for the possibility that the PDB file has atoms in a non-standard order, or has different subsets of atoms, such as all atoms including hydrogens, while MODELLER uses only heavy atoms, or vice versa). We then create a selection of all atoms, since most MODELLER energy functions can operate on a subset of model atoms. The DOPE energy is then calculated with the assess_dope command, and we additionally request an energy profile, smoothed over a 15 residue window, and normalized by the number of restraints acting on each residue. This profile is written to a file "TvLDH.profile", which can be used as input to a graphing program. For example, it could be plotted with GNUPLOT using the command 'plot "TvLDH.profile" using 1:42 with lines'. Alternatively, you can use the plot_profiles.py script included in the tutorial zip file to plot profiles with the Python matplotlib package.

The GA341 score, as well as external analysis with the PROCHECK program, confirms that TvLDH.B99990001.pdb is a reasonable model. However, the plotted DOPE score profile (below) shows regions of relatively high energy for the long active site loop between residues 90 and 100 and the long helices at the C-terminal end of the target sequence. (Note that we have superposed the model profile on the template profile - gaps in the plot can be seen corresponding to the gaps in the alignment. Remember that the scores are not absolute, so we cannot make a direct numerical comparison between the two. However, we can get an idea of the quality of our input alignment this way by comparing the rough shapes of the two profiles - if one is obviously shifted relative to the other, it is likely that the alignment is also shifted from the correct one.)




# Advanced Modeling. Model a sequence based on multiple templates and bound to a ligand
The structure of the Malate Dehydrogenase 1bdm has been clustered in the DBAli database within the family fm00495 of 4 members (2mdh:A, 2mdh:B. 1b8p:A and 1bdm:A). The multiple alignment generated by the command salign() in MODELLER is used in DBAli to generate a multiple structure alignment of the family. The alignment can be downloaded from the DBAli database or you can use the file `salign.py' to calculate it on your computer.

# Illustrates the SALIGN multiple structure/sequence alignment

from modeller import *

log.verbose()
env = environ()
env.io.atom_files_directory = './:../atom_files/'

aln = alignment(env)
for (code, chain) in (('2mdh', 'A'), ('1bdm', 'A'), ('1b8p', 'A')):
    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               dendrogram_file='fm00495.tree',
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

aln.write(file='fm00495.pap', alignment_format='PAP')
aln.write(file='fm00495.ali', alignment_format='PIR')

aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
           rr_file='$(LIB)/as1.sim.mat', overhang=30,
           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
           gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
           alignment_type='progressive', feature_weights=[0]*6,
           improve_alignment=False, fit=False, write_fit=True,
           write_whole_pdb=False, output='QUALITY')

File: multiple_template/salign.py

The reads in all of the sequences from PDB files (using the append_model command), and then uses salign multiple times, to generate an initial rough alignment and then improve upon it by using more information. The alignment is then written out in both PIR and PAP formats, and a quality score is calculated by calling salign one more time.

After inspecting the multiple structure alignment it is evident that chain B of 2mdh contains an unusual number of LYS residues. The HEADER of the PDB file indicates that the sequence of the protein was unknown at the time of refinement and it was difficult to identify most of the residues in the structure. Therefore, the 2mdh:B entry was removed from the multiple structure alignment.

 _aln.pos         10        20        30        40        50        60
2mdhA     GSMQIRVLVTG-AAQLAFTLLYSIGDGSVFGKNQPILLSLMDVVP--KQQTSEAVNMQLQNCALP-LL 
1bdmA     MKAPVRVAVTGAAGQIGYSLLFRIAAGEMLGKDQPVILQLLEIPQ--AMKALEGVVMELEDCAFPLLA 
1b8pA     -KTPMRVAVTGAAGQICYSLLFRIANGDMLGKDQPVILQLLEIPNEKAQKALQGVMMEIDDCAFPLLA 
 _consrvd      ** *** * *    **  *  *   ** **  * *              * *    ** * *

 _aln.p   70        80        90       100       110       120       130
2mdhA     KSQFGKNSGN-YASQNVGVLLAGQRAKNAAKN---LKANVKIFKCQGAALNKYWKKSVIVIVVGNPAT 
1bdmA     GLEATDDPDVAFKDADYALLVGAAPR---------LQVNGKIFTEQGRALAEVAKKDVKVLVVGNPAN 
1b8pA     GMTAHADPMTAFKDADVALLVGARPRGPGMERKDLLEANAQIFTVQGKAIDAVASRNIKVLVVGNPAN 
 _consrvd                    *               *  *  **  ** *          * ******

 _aln.pos  140       150       160       170       180       190       200
2mdhA     NNCLTASKNSAQLNKAKQVNSVKLNHNRAKSMLSQKLGNSPKLSKNVILYGQHGQSQFSGLIQLQLQN 
1bdmA     TNALIAYKNAPGLNPRNFTAMTRLDHNRAKAQLAKKTGTGVDRIRRMTVWGNHSSIMFPDLFHAEVD- 
1b8pA     TNAYIAMKSAPSLPAKNFTAMLRLDHNRALSQIAAKTGKPVSSIEKLFVWGNHSPTMYADYRYAQID- 
 _consrvd  *   * *    *          * ****      * *            * *

 _aln.pos    210       220       230       240       250       260       270
2mdhA     KQSAGVR-ASKNQSWKTSIYNNVIQQRGVVHVQARTANNSMKTGFALNLYVKHLWKGISQ-KLAQMGL 
1bdmA     --GRPALELV-DMEWYEKVFIPTVAQRGAAIIQARGASSAASAANAAIEHIRDWALGTPEGDWVSMAV 
1b8pA     --GASVKDMINDDAWNRDTFLPTVGKRGAAIIDARGVSSAASAANAAIDHIHDWVLGTAG-KWTTMGI 
 _consrvd               *           **     **          *          *        *

 _aln.pos      280       290       300       310       320       330
2mdhA     IAHGKAAASPKQNFSCVTRLQNKTWKIVEGLPINDFSREKMNETAKELAEEETEFAEKNSNA 
1bdmA     PSQGEYGIPEGIVYSFPVTAKDGAYRVVEGLEINEFARKRMEITAQELLDEMEQVKALGLI- 
1b8pA     PSDGSYGIPEGVIFGFPVTTENGEYKIVQGLSIDAFSQERINVTLNELLEEQNGVQ-HLLG- 
 _consrvd    *                       * ** *  *       *  **  *
File: multiple_template/fm00495.pap

As for the basic example in the tutorial, next we need to align our query sequence to the template structures. For that task we again use the salign() command (file `align2d_mult.py'). We set the align_block parameter to equal the number of structures in the template alignment, len(aln), (i.e. 3), and request a pairwise alignment, since we do not want to change the existing alignment between the templates. By setting gap_function we request the use of a structure-dependent gap penalty, using structural information for these 3 sequences. Only sequence information is used for the final TvLDH sequence.

from modeller import *

log.verbose()
env = environ()

env.libs.topology.read(file='$(LIB)/top_heav.lib')

# Read aligned structure(s):
aln = alignment(env)
aln.append(file='fm00495.ali', align_codes='all')
aln_block = len(aln)

# Read aligned sequence(s):
aln.append(file='TvLDH.ali', align_codes='TvLDH')

# Structure sensitive variable gap penalty sequence-sequence alignment:
aln.salign(output='', max_gap_length=20,
           gap_function=True,   # to use structure-dependent gap penalty
           alignment_type='PAIRWISE', align_block=aln_block,
           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
           gap_penalties_1d=(-450, 0),
           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
           similarity_flag=True)

aln.write(file='TvLDH-mult.ali', alignment_format='PIR')
aln.write(file='TvLDH-mult.pap', alignment_format='PAP')
File: multiple_template/align2d_mult.py

Next, we build the new model for the TvLDH target sequence based on the alignment against the multiple templates using the `model_mult.py' file:

from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='TvLDH-mult.ali',
              knowns=('1bdmA','2mdhA','1b8pA'), sequence='TvLDH')
a.starting_model = 1
a.ending_model = 5
a.make()
File: multiple_template/model_mult.py

Finally, we use the DOPE potential to evaluate the new model coordinates using the `evaluate_model.py' file:

from modeller import *
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# read model file
mdl = complete_pdb(env, 'TvLDH.B99990001.pdb')

# Assess all atoms with DOPE:
s = selection(mdl)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='TvLDH.profile',
              normalize_profile=True, smoothing_window=15)
File: multiple_template/evaluate_model.py

The evaluation of the model indicates that the problematic loop (residues 90 to 100) has improved by using multiple structural templates. The global DOPE score for the models also improved from -38999.7 to -39164.4. MODELLER was able to use the variability in the loop region from the three templates to generate a more accurate conformation of the loop. However, the conformation of a loop in the region around the residue 275 at the C-terminal end of the sequence has higher DOPE score than for the model based on a single template.


We will use the loopmodel class in MODELLER to refine the conformation of the loop between residues 273 and 283. We will use the model number 1 created in the previous example as a starting structure to refine the loop. You can find this structure renamed as `TvDLH_mult.pdb' in the loop_modeling subdirectory.

Loop refining

The loop optimization method relies on a scoring function and optimization schedule adapted for loop modeling. It is used automatically to refine comparative models if you use the loopmodel class rather than automodel; see the example below.
# Loop refinement of an existing model
from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = './:../atom_files'

# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # 10 residue insertion 
        return selection(self.residue_range('273', '283'))

m = MyLoop(env,
           inimodel='TvLDH-mult.pdb', # initial model of the target
           sequence='TvLDH')          # code of the target

m.loop.starting_model= 1           # index of the first loop model 
m.loop.ending_model  = 10          # index of the last loop model
m.loop.md_level = refine.very_fast # loop refinement method; this yields
                                   # models quickly but of low quality;
                                   # use refine.slow for better models

m.make()

File: loop_modeling/loop_refine.py

In this example, the loopmodel class is used to refine a region of an existing coordinate file. Note that this example also redefines the loopmodel.select_loop_atoms routine. This is necessary in this case, as the default selection selects all gaps in the alignment for refinement, and in this case no alignment is available. You can still redefine the routine, even if you do have an alignment, if you want to select a different region for optimization. Note that for the sake of time, we will be building only 10 different independently optimized loop conformations by setting the loop.ending_model parameter to 10. The next image shows the superimposition of the 10 conformations of the loop modeling. In blue, green and red we have marked the initial, best and worst loop conformations as scored by DOPE, respectively.


The file `model_energies.py' computes the DOPE score for all built models by using a Python for loop. The best energy loop corresponds to the 8th model (file: `model_energies.py') with a global DOPE score of -39099.1. Its energy profile calculated by `evaluate_model.py' is shown next.


here is only a very small increase of global DOPE score by ab-initio refinement of the loop. However, there is a small decrease in the DOPE score in the region of the loop. Therefore, we will continue the next step using the best refined structure (file: `TvLDH.BL00080001.pdb'), which is renamed in the ligand directory as `TvLDH-loop.pdb'. It is important to note that a most accurate approach to loop refinement requires the modeling of hundreds of independent conformations and their clustering to select the most representative structures of the loop.

Modeling ligands in the binding site

1emd, a malate dehydrogenase from E. coli, was identified in PDB. While the 1emd sequence shares only 32% sequence identity with TvLDH, the active site loop and its environment are more conserved. The loop for residues 90 to 100 in the 1emd structure is well resolved. Moreover, 1emd was solved in the presence of a citrate substrate analog and the NADH cofactor. The new alignment in the PAP format is shown below (file `TvLDH-1emd_bs.pap').

 _aln.pos            10        20        30        40        50        60
TvLDH        MSEAAHVLITGAAGQIGYILSHWIASGELYGDRQVYLHLLDIPPAMNRLTALTMELEDCAFPHLA 
TvLDH_model  MSEAAHVLITGAAGQIGYILSHWIASGELYGDRQVYLHLLDIPPAMNRLTALTMELEDCAFPHLA 
1emd         ----------------------------------------------------------------- 
 _consrvd


 _aln.pos       70        80        90       100       110       120       130
TvLDH        GFVATTDPKAAFKDIDCAFLVASMPLKPGQVRADLISSNSVIFKNTGEYLSKWAKPSVKVLVIGN 
TvLDH_model  GFVATTDPKAAFKDIDCAFLVASMPLKPGQVRADLISSNSVIFKNTGEYLSKWAKPSVKVLVIGN 
1emd         ----------------------GVRRKPGMDRSDLFNVN-------------------------- 
 _consrvd                              ***  * **   *


 _aln.pos           140       150       160       170       180       190
TvLDH        PDNTNCEIAMLHAKNLKPENFSSLSMLDQNRAYYEVASKLGVDVKDVHDIIVWGNHGESMVADLT 
TvLDH_model  PDNTNCEIAMLHAKNLKPENFSSLSMLDQNRAYYEVASKLGVDVKDVHDIIVWGNHGESMVADLT 
1emd         ----------------------------------------------------------------- 
 _consrvd


 _aln.pos      200       210       220       230       240       250       260
TvLDH        QATFTKEGKTQKVVDVLDHDYVFDTFFKKIGHRAWDILEHRGFTSAASPTKAAIQHMKAWLFGTA 
TvLDH_model  QATFTKEGKTQKVVDVLDHDYVFDTFFKKIGHRAWDILEHRGFTSAASPTKAAIQHMKAWLFGTA 
1emd         ----------------------------------------------------------------- 
 _consrvd


 _aln.pos           270       280       290       300       310       320
TvLDH        PGEVLSMGIPVPEGNPYGIKPGVVFSFPCNVDKEGKIHVVEGFKVNDWLREKLDFTEKDLFHEKE 
TvLDH_model  PGEVLSMGIPVPEGNPYGIKPGVVFSFPCNVDKEGKIHVVEGFKVNDWLREKLDFTEKDLFHEKE 
1emd         ----------------------------------------------------------------- 
 _consrvd


 _aln.pos      330
TvLDH        IALNHLAQGG/.. 
TvLDH_model  IALNHLAQGG/-- 
1emd         ----------/.. 
 _consrvd
File: ligand/TvLDH-1emd_bs.pap

The modified alignment refers to an edited 1emd structure (1emd_bs), as a second template. The alignment corresponds to a model that is based on 1emd_bs in its active site loop and on TvLDH_model, which corresponds to the best model from the previous step, in the rest of the fold. Four residues on both sides of the active site loop are aligned with both templates to ensure that the loop has a good orientation relative to the rest of the model.

The modeling script below has several changes with respect to `model-single.py'. First, the name of the alignment file assigned to alnfile is updated. Next, the variable knowns is redefined to include both templates. Another change is an addition of the `env.io.hetatm = True' command to allow reading of the non-standard pyruvate and NADH residues from the input PDB files. The script is shown next (file `model-multiple-hetero.py').

from modeller import *
from modeller.automodel import *

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        for ids in (('NH1:161:A', 'O1A:336:B'),
                    ('NH2:161:A', 'O1B:336:B'),
                    ('NE2:186:A', 'O2:336:B')):
            atoms = [self.atoms[i] for i in ids]
            rsr.add(forms.upper_bound(group=physical.upper_distance,
                                      feature=features.distance(*atoms),
                                      mean=3.5, stdev=0.1))

env = environ()
env.io.hetatm = True
a = MyModel(env, alnfile='TvLDH-1emd_bs.ali',
              knowns=('TvLDH_model','1emd'), sequence='TvLDH')
a.starting_model = 1
a.ending_model = 5
a.make()
File: ligand/model-multiple-hetero.py

A ligand can be included in a model in two ways by MODELLER. The first case corresponds to the ligand that is not present in the template structure, but is defined in the MODELLER residue topology library. Such ligands include water molecules, metal ions, nucleotides, heme groups, and many other ligands (see question 8 in the the MODELLER FAQ). This situation is not explored further here. The second case corresponds to the ligand that is already present in the template structure. We can assume either that the ligand interacts similarly with the target and the template, in which case we can rely on MODELLER to extract and satisfy distance restraints automatically, or that the relative orientation is not necessarily conserved, in which case the user needs to supply restraints on the relative orientation of the ligand and the target (the conformation of the ligand is assumed to be rigid). The two cases are illustrated by the NADH cofactor and pyruvate modeling, respectively. Both NADH and cofactor are indicated by the `.' characters at the end of each sequence in the alignment file above (the `/' character indicates a chain break). In general, the `.' character in MODELLER indicates an arbitrary generic residue called a ``block'' residue (for details see the section on block residues in the MODELLER manual). Note that the `.' characters are present both in one of the template structures and in the model sequence. The former tells MODELLER to read the ligands from the template, and the latter tells it to include the ligands in the model. The 1emd structure file contains a citrate substrate analog. To obtain a model with pyruvate, the physiological substrate of TvLDH, we convert the citrate analog in 1emd into pyruvate by deleting the group CH(COOH)2, thus obtaining the 1emd_bs template file. A major advantage of using the `.' characters is that it is not necessary to define the residue topology.

To obtain the restraints on pyruvate, we first superpose the structures of several LDH and MDH enzymes solved with ligands. Such a comparison allows us to identify absolutely conserved electrostatic interactions involving catalytic residues Arg161 and His186 on one hand, and the oxo groups of the lactate and malate ligands on the other hand. The modeling script can now be expanded by creating a new class 'MyModel', which is derived from automodel but which differs in one important respect: the special_restraints routine is redefined to add, to the default restraints, some user defined distance restraints between the conserved atoms of the active site residues and their substrate. In this case, a harmonic upper bound restraint of 3.5±0.1Å is imposed on the distances between the three specified pairs of atoms. A trick is used to prevent MODELLER from automatically calculating distance restraints on the pyruvate-TvLDH complex; the ligand in the 1emd_bs template is moved beyond the upper bound on the ligand-protein distance restraints (i.e., 10).

The final selected model (shown in the ribbons image below) has a DOPE global score of -37640.9. The DOPE score is increased due to the new interactions of the protein with the ligand that are not accounted when calculating the DOPE score.


#Iterative Modeling. Increase the accuracy of the modeling exercise by iterating the 4 step process.

Several structures of dihydrofolate reductase (DHFR) are known. However, the structure of DHFR from Haloferax volcanii was not known and its sequence identity with DHFRs of known structure is rather low at ~30%. A model of H. volcanii DHFR (HVDFR) was constructed before the experimental structure was solved. This example illustrates the power of the iterative alignment-modeling-evaluation approach to comparative modeling.

Of all the available DHFR structures, HVDHFR has the sequence most similar to DHFR from E. coli. The PDB entry 4dfr corresponds to a high resolution (1.7Å) E. coli DHFR structure. It contains two copies of the molecule, named chain A and chain B. According to the authors, the structure for chain B is of better quality than that of chain A. The following script file aligns HVDFR and chain B of 4dfr.

from modeller import *

env = environ()
mdl = model(env, file='4dfr.pdb', model_segment=('FIRST:B', 'LAST:B'))
aln = alignment(env)
aln.append_model(mdl, align_codes='4dfr')
aln.append(file='hvdfr.seq', align_codes='hvdfr')
aln.align2d()
aln.write(file='hvdfr-4dfr.ali')
aln.write(file='hvdfr-4dfr.pap', alignment_format='PAP',
          alignment_features='INDICES HELIX BETA')
File: align2d-4.py

Some options used in this example include model_segment, which is used to indicate chain B of 4dfr; and alignment_features, which is used to output information such as secondary structure to the alignment file in the PAP format.

 _aln.pos         10        20        30        40        50        60
4dfr      -MISLIAALAVDRVIGMENAMPW-NLPADLAWFKRNTLDKPVIMGRHTWESIGRPLPGRKNIILSSQP 
hvdfr     MELVSVAALAENRVIGRDGELPWPSIPADKKQYRSRIADDPVVLGRTTFESMRDDLPGSAQIVMSRSE 
 _helix                            999999999999       999999999
 _beta     9999999999                           999999             99999999


 _aln.p   70        80        90       100       110       120       130
4dfr      GTDDRVTWVKSV----DEAIAACGDVPEIMVIGGGRVYEQFLPKAQKLYLTHIDAEVEGDTHFPDYEP 
hvdfr     RSFSVDTAHRAASVEEAVDIAASLDAETAYVIGGAAIYALFQPHLDRMVLSRVPGEYEGDTYYPEWDA 
 _helix             99    99999999         99999999
 _beta         99999                9999999            999999999


 _aln.pos  140       150       160
4dfr      DDWESVFSEFHDADAQNSHSYCFKILERR 
hvdfr     AEWELDAETDHEG---FTLQEWVRSASSR 
 _helix
 _beta     999999999999    999999999999
File: hvdfr-4dfr.pap

Using the PIR alignment file hvdfr-4dfr.ali, an initial model is calculated:

from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='hvdfr-4dfr.ali',
              knowns='4dfr', sequence='hvdfr')
a.starting_model = 1
a.ending_model = 1
a.make()
File: model4.py

Because the sequence identity between 4dfr and HVDFR is relatively low (30%), the automated alignment is likely to contain errors. The evaluation of the model with the DOPE potential in MODELLER shows two regions with positive energy.

PROSAII profile for model initial model

DOPE score profile for model inihvdfr.B99990001.pdb

The first region is around residue 85, the second region is at the C-terminal end of the protein. Referring to the target--template alignment above (hvdfr-4dfr.pap), it is easy to understand why the first positive peak appears. The insertion around position 85 of the alignment was placed in the middle of an α-helix in the template (the "9" characters on the first line below the sequence mark the helices). Moving the insertion to the end of the α-helix may improve the model.

The second problem, which occurs in the C-terminal region of the alignment, is less clear. The deletion in that region of the alignment corresponds to the loop between the last two β-strands of 4dfr (a β-hairpin). Since the profile suggests that this region is in error, an alternative alignment should be tried. One possibility is that the deletion is actually longer, making the C-terminal β-hairpin shorter in HVDFR. One plausible alignment based on these considerations is shown here.

 _aln.pos         10        20        30        40        50        60
4dfr      M-ISLIAALAVDRVIGMENAMPW-NLPADLAWFKRNTLDKPVIMGRHTWESIGRPLPGRK 
hvdfr     MELVSVAALAENRVIGRDGELPWPSIPADKKQYRSRIADDPVVLGRTTFESMRDDLPGSA 
 _helix                            999999999999       999999999
 _beta    9 999999999                           999999             999


 _aln.pos         70        80        90       100       110       120
4dfr      NIILSSQPGT--DDRVTWVKSVDEAIAACG--DVPEIMVIGGGRVYEQFLPKAQKLYLTH 
hvdfr     QIVMSRSERSFSVDTAHRAASVEEAVDIAASLDAETAYVIGGAAIYALFQPHLDRMVLSR 
 _helix                       9999999999           99999999
 _beta    99999          99999              9999999            9999999


 _aln.pos        130       140       150       160
4dfr      IDAEVEGDTHFPDYEPDDWESVFSEFHDADAQNSHSYCFKILERR---- 
hvdfr     VPGEYEGDTYYPEWDAAEWELDAETDHE-------GFTLQEWVRSASSR 
 _helix
 _beta    99               999999999999    999999999999
File: hvdfr-4dfr-2.pap

A new model was calculated using this alignment and the script, modified to use the new alignment (see file `model5.py'). Its DOPE score profile is shown in the next figure.

DOPE profile for model hvdfr.B99990001.pdb

DOPE profile for model hvdfr.B99990001.pdb

The main positive peak disappeared and the new global DOPE score dropped from -15498.7 to -15720.3.

The process outlined here could be iterated until no improvement in the evaluation can be achieved. This iterative alignment-building-evaluation process has been developed in the MOULDER protocol which will be further implemented in MODELLER .

# Difficult Modeling. Model a sequence based on a low identity to a template

Several structures of dihydrofolate reductase (DHFR) are known. However, the structure of DHFR from Haloferax volcanii was not known and its sequence identity with DHFRs of known structure is rather low at ~30%. A model of H. volcanii DHFR (HVDFR) was constructed before the experimental structure was solved. This example illustrates the power of the iterative alignment-modeling-evaluation approach to comparative modeling.

Of all the available DHFR structures, HVDHFR has the sequence most similar to DHFR from E. coli. The PDB entry 4dfr corresponds to a high resolution (1.7Å) E. coli DHFR structure. It contains two copies of the molecule, named chain A and chain B. According to the authors, the structure for chain B is of better quality than that of chain A. The following script file aligns HVDFR and chain B of 4dfr.

from modeller import *

env = environ()
mdl = model(env, file='4dfr.pdb', model_segment=('FIRST:B', 'LAST:B'))
aln = alignment(env)
aln.append_model(mdl, align_codes='4dfr')
aln.append(file='hvdfr.seq', align_codes='hvdfr')
aln.align2d()
aln.write(file='hvdfr-4dfr.ali')
aln.write(file='hvdfr-4dfr.pap', alignment_format='PAP',
          alignment_features='INDICES HELIX BETA')
File: align2d-4.py

Some options used in this example include model_segment, which is used to indicate chain B of 4dfr; and alignment_features, which is used to output information such as secondary structure to the alignment file in the PAP format.

 _aln.pos         10        20        30        40        50        60
4dfr      -MISLIAALAVDRVIGMENAMPW-NLPADLAWFKRNTLDKPVIMGRHTWESIGRPLPGRKNIILSSQP 
hvdfr     MELVSVAALAENRVIGRDGELPWPSIPADKKQYRSRIADDPVVLGRTTFESMRDDLPGSAQIVMSRSE 
 _helix                            999999999999       999999999
 _beta     9999999999                           999999             99999999


 _aln.p   70        80        90       100       110       120       130
4dfr      GTDDRVTWVKSV----DEAIAACGDVPEIMVIGGGRVYEQFLPKAQKLYLTHIDAEVEGDTHFPDYEP 
hvdfr     RSFSVDTAHRAASVEEAVDIAASLDAETAYVIGGAAIYALFQPHLDRMVLSRVPGEYEGDTYYPEWDA 
 _helix             99    99999999         99999999
 _beta         99999                9999999            999999999


 _aln.pos  140       150       160
4dfr      DDWESVFSEFHDADAQNSHSYCFKILERR 
hvdfr     AEWELDAETDHEG---FTLQEWVRSASSR 
 _helix
 _beta     999999999999    999999999999
File: hvdfr-4dfr.pap

Using the PIR alignment file hvdfr-4dfr.ali, an initial model is calculated:

from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='hvdfr-4dfr.ali',
              knowns='4dfr', sequence='hvdfr')
a.starting_model = 1
a.ending_model = 1
a.make()
File: model4.py

Because the sequence identity between 4dfr and HVDFR is relatively low (30%), the automated alignment is likely to contain errors. The evaluation of the model with the DOPE potential in MODELLER shows two regions with positive energy.

PROSAII profile for model initial model

DOPE score profile for model inihvdfr.B99990001.pdb

The first region is around residue 85, the second region is at the C-terminal end of the protein. Referring to the target--template alignment above (hvdfr-4dfr.pap), it is easy to understand why the first positive peak appears. The insertion around position 85 of the alignment was placed in the middle of an α-helix in the template (the "9" characters on the first line below the sequence mark the helices). Moving the insertion to the end of the α-helix may improve the model.

The second problem, which occurs in the C-terminal region of the alignment, is less clear. The deletion in that region of the alignment corresponds to the loop between the last two β-strands of 4dfr (a β-hairpin). Since the profile suggests that this region is in error, an alternative alignment should be tried. One possibility is that the deletion is actually longer, making the C-terminal β-hairpin shorter in HVDFR. One plausible alignment based on these considerations is shown here.

 _aln.pos         10        20        30        40        50        60
4dfr      M-ISLIAALAVDRVIGMENAMPW-NLPADLAWFKRNTLDKPVIMGRHTWESIGRPLPGRK 
hvdfr     MELVSVAALAENRVIGRDGELPWPSIPADKKQYRSRIADDPVVLGRTTFESMRDDLPGSA 
 _helix                            999999999999       999999999
 _beta    9 999999999                           999999             999


 _aln.pos         70        80        90       100       110       120
4dfr      NIILSSQPGT--DDRVTWVKSVDEAIAACG--DVPEIMVIGGGRVYEQFLPKAQKLYLTH 
hvdfr     QIVMSRSERSFSVDTAHRAASVEEAVDIAASLDAETAYVIGGAAIYALFQPHLDRMVLSR 
 _helix                       9999999999           99999999
 _beta    99999          99999              9999999            9999999


 _aln.pos        130       140       150       160
4dfr      IDAEVEGDTHFPDYEPDDWESVFSEFHDADAQNSHSYCFKILERR---- 
hvdfr     VPGEYEGDTYYPEWDAAEWELDAETDHE-------GFTLQEWVRSASSR 
 _helix
 _beta    99               999999999999    999999999999
File: hvdfr-4dfr-2.pap

A new model was calculated using this alignment and the script, modified to use the new alignment (see file `model5.py'). Its DOPE score profile is shown in the next figure.

DOPE profile for model hvdfr.B99990001.pdb

DOPE profile for model hvdfr.B99990001.pdb

The main positive peak disappeared and the new global DOPE score dropped from -15498.7 to -15720.3.

The process outlined here could be iterated until no improvement in the evaluation can be achieved. This iterative alignment-building-evaluation process has been developed in the MOULDER protocol which will be further implemented in MODELLER .


# Modeling with cryo-EM. Model a sequence using both template and cryo-EM data.
The use of MODELLER to build complete 3D models of an 'unknown' protein structure (i.e. not solved by X-ray crystallography) given three sources of experimental information:

The amino acid sequence of the unknown protein structure.
The PDB database of other, known protein structures.
A cryo-EM map of the unknown protein structure.















