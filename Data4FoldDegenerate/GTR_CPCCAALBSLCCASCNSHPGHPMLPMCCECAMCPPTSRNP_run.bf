/* T. Latrille.
Hyphy inference for an "experimental" dataset
*/
global mu=1; global pnC=0.25; global pnG=0.25; global pnA=0.25; global exchCT=1.0; global exchCG=1.0; global exchAG=1.0; global exchGT=1.0; global exchAT=1.0; 
global pnT:=1.0-(pnC+pnG+pnA); 

RANDOM_STARTING_PERTURBATIONS = 1;

StateMatrix={4, 4, 
{0,1,mu*t*pnC} {0,2,mu*t*exchAG*pnG} {0,3,mu*t*exchAT*pnT} {1,0,mu*t*pnA} {1,2,mu*t*exchCG*pnG} {1,3,mu*t*exchCT*pnT} {2,0,mu*t*exchAG*pnA} {2,1,mu*t*exchCG*pnC} {2,3,mu*t*exchGT*pnT} {3,0,mu*t*exchAT*pnA} {3,1,mu*t*exchCT*pnC} {3,2,mu*t*exchGT*pnG} }
;

/* Read in the data */
DataSet	raw_data = ReadDataFile("/home/thibault/NucleotideBias/Data4FoldDegenerate/CDS.CPCCAALBSLCCASCNSHPGHPMLPMCCECAMCPPTSRNP.fasta");


DataSetFilter filt_data = CreateFilter(raw_data,1);


/* Set up frequencies.  */
StateFrequencies={{pnA},{pnC},{pnG},{pnT}};;

/* Optimize likelihoods for each frequency specification */
////////////// r5_f3_w0 //////////////
global mu=1; global pnC=0.25; global pnG=0.25; global pnA=0.25; global exchCT=1.0; global exchCG=1.0; global exchAG=1.0; global exchGT=1.0; global exchAT=1.0; 
global pnT:=1.0-(pnC+pnG+pnA); 

Model StateModel = (StateMatrix, StateFrequencies, 0);
UseModel (USE_NO_MODEL);
UseModel(StateModel);
Tree    Tree_newick = (((Callicebus_donacophilus:0.0147108,(Pithecia_pithecia:0.00952312,(Chiropotes_chiropotes:0.00571323,Cacajao_calvus:0.00533583)ChiropoCacajao:0.003026)PitheChiroCacaj:0.00472266)CallPithChirCaca:0.00223835,((Alouatta_palliata:0.0101252,(Ateles_belzebuth:0.00650514,(Lagothrix_cana:0.00434368,Brachyteles_arachnoides:0.00767233)LagothrBrachyt:0.00066753)AteleLagotBrach:0.00186522)AlouAtelLagoBrac:0.00456489,((Saguinus_fuscicollis:0.00928601,(Leontopithecus_rosalia:0.0126459,(Callimico_goeldii:0.00772146,Callithrix_jacchus:0.0081631)CallimiCallith:0.00279215)LeontCalliCalli:0.00081491)SaguLeonCallCall:0.0061857,(Aotus_azarai:0.0113188,(Saimiri_oerstedii:0.012929,Cebus_apella:0.0106798)SaimiriCebus_a:0.00191902)AotusSaimiCebus:0.00037732)SaLeCaCaAoSaCe:0.00267581)AlAtLaBrSaLeCaCaAoSaCe:0.00094259)CPCCAALBSLCCASC:0.0184643,(((Nomascus_siki:0.00235839,(Symphalangus_syndactylus:0.00282525,Hylobates_agilis:0.00338925)SymphalHylobat:0.00031551)NomasSymphHylob:0.00627798,(Pongo_pygmaeus:0.00765437,(Gorilla_gorilla:0.00332432,(Homo_sapiens:0.00287512,Pan_paniscus:0.00251746)Homo_saPan_pan:0.0005525)GorilHomo_Pan_p:0.00362861)PongGoriHomoPan_:0.00144949)NoSyHyPoGoHoPa:0.00377405,(((Macaca_mulatta:0.00318493,((Lophocebus_aterrimus:0.00310964,Papio_papio:0.0012536)LophocePapio_p:0.00135085,(Mandrillus_sphinx:0.00235559,Cercocebus_torquatus:0.00150973)MandrilCercoce:0.0005635)LophPapiMandCerc:0.00041727)MacLopPapManCer:0.00177008,((Chlorocebus_aethiops:0.00188392,Erythrocebus_patas:0.0029496)ChlorocErythro:0.00105524,(Cercopithecus_albogularis:0.00208713,(Allenopithecus_nigroviridis:0.00348117,Miopithecus_ogouensis:0.00285094)AllenopMiopith:0.00094893)CercoAllenMiopi:5.238e-05)ChlEryCerAllMio:0.00099544)MaLoPaMaCeChErCeAlMi:0.00280935,((Colobus_angolensis:0.00458254,Piliocolobus_badius:0.0044103)ColobusPilioco:0.00122269,(Presbytis_melalophos:0.00297883,((Trachypithecus_francoisi:0.00079858,Semnopithecus_entellus:0.00163174)TrachypSemnopi:0.00264889,(Rhinopithecus_brelichi:0.00262893,(Nasalis_larvatus:0.0022615,Pygathrix_cinerea:0.00273209)NasalisPygathr:0.0001936)RhinoNasalPygat:0.00068921)TraSemRhiNasPyg:0.00012943)PreTraSemRhiNasPyg:0.00102253)CoPiPrTrSeRhNaPy:0.00227456)MLPMCCECAMCPPTSRNP:0.00852039)NSHPGHPMLPMCCECAMCPPTSRNP:0.00738331);;

branchNames 	= BranchName (Tree_newick, -1);
branchLengths	= BranchLength (Tree_newick, -1);

for (k = 0; k < Columns(branchNames)-1; k=k+1)
{
	ExecuteCommands("Tree_newick." + branchNames[k] + ".t:=" + branchLengths[k] + ";");
}

LikelihoodFunction  LikStateModel = (filt_data, Tree_newick);
Optimize (paramValues, LikStateModel);
fprintf (stdout, LikStateModel);
fprintf ("GTR_CPCCAALBSLCCASCNSHPGHPMLPMCCECAMCPPTSRNP_run.bf_hyout.txt", LikStateModel);