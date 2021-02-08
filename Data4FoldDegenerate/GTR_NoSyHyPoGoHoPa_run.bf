/* T. Latrille.
Hyphy inference for an "experimental" dataset
*/
global mu=1; global pnG=0.25; global pnA=0.25; global pnC=0.25; global exchCT=1.0; global exchAG=1.0; global exchAT=1.0; global exchGT=1.0; global exchCG=1.0; 
global pnT:=1.0-(pnG+pnA+pnC); 

RANDOM_STARTING_PERTURBATIONS = 1;

StateMatrix={4, 4, 
{0,1,mu*t*pnC} {0,2,mu*t*exchAG*pnG} {0,3,mu*t*exchAT*pnT} {1,0,mu*t*pnA} {1,2,mu*t*exchCG*pnG} {1,3,mu*t*exchCT*pnT} {2,0,mu*t*exchAG*pnA} {2,1,mu*t*exchCG*pnC} {2,3,mu*t*exchGT*pnT} {3,0,mu*t*exchAT*pnA} {3,1,mu*t*exchCT*pnC} {3,2,mu*t*exchGT*pnG} }
;

/* Read in the data */
DataSet	raw_data = ReadDataFile("/home/thibault/NucleotideBias/Data4FoldDegenerate/CDS.NoSyHyPoGoHoPa.fasta");


DataSetFilter filt_data = CreateFilter(raw_data,1);


/* Set up frequencies.  */
StateFrequencies={{pnA},{pnC},{pnG},{pnT}};;

/* Optimize likelihoods for each frequency specification */
////////////// r5_f3_w0 //////////////
global mu=1; global pnG=0.25; global pnA=0.25; global pnC=0.25; global exchCT=1.0; global exchAG=1.0; global exchAT=1.0; global exchGT=1.0; global exchCG=1.0; 
global pnT:=1.0-(pnG+pnA+pnC); 

Model StateModel = (StateMatrix, StateFrequencies, 0);
UseModel (USE_NO_MODEL);
UseModel(StateModel);
Tree    Tree_newick = ((Nomascus_siki:0.00235839,(Symphalangus_syndactylus:0.00282525,Hylobates_agilis:0.00338925)SymphalHylobat:0.00031551)NomasSymphHylob:0.00627798,(Pongo_pygmaeus:0.00765437,(Gorilla_gorilla:0.00332432,(Homo_sapiens:0.00287512,Pan_paniscus:0.00251746)Homo_saPan_pan:0.0005525)GorilHomo_Pan_p:0.00362861)PongGoriHomoPan_:0.00144949);;

branchNames 	= BranchName (Tree_newick, -1);
branchLengths	= BranchLength (Tree_newick, -1);

for (k = 0; k < Columns(branchNames)-1; k=k+1)
{
	ExecuteCommands("Tree_newick." + branchNames[k] + ".t:=" + branchLengths[k] + ";");
}

LikelihoodFunction  LikStateModel = (filt_data, Tree_newick);
Optimize (paramValues, LikStateModel);
fprintf (stdout, LikStateModel);
fprintf ("GTR_NoSyHyPoGoHoPa_run.bf_hyout.txt", LikStateModel);