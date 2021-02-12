/* T. Latrille.
Hyphy inference for an "experimental" dataset
*/
{% for var in vars_list %}{{ var }} {% endfor %}
{% for var in constrains_list %}{{ var }} {% endfor %}

StateMatrix={{ matrix }};

/* Read in the data */
DataSet	raw_data = ReadDataFile("{{ fasta_infile }}");

{% if filter_codons_stop %}
/* Filter the data to find and remove any stop codons*/
DataSetFilter filt_data = CreateFilter(raw_data,3,"","","TAA,TAG,TGA");
{% else %}
DataSetFilter filt_data = CreateFilter(raw_data,1);
{% endif %}

/* Set up frequencies.  */
StateFrequencies={{ freqs }};

/* Optimize likelihoods for each frequency specification */
////////////// {{ param }} //////////////
{% for var in vars_list %}{{ var }} {% endfor %}
{% for var in constrains_list %}{{ var }} {% endfor %}

Model StateModel = (StateMatrix, StateFrequencies, 0);
UseModel (USE_NO_MODEL);
UseModel(StateModel);
Tree    Tree_newick = {{ tree }};

branchNames 	= BranchName (Tree_newick, -1);
branchLengths	= BranchLength (Tree_newick, -1);

for (k = 0; k < Columns(branchNames)-1; k=k+1)
{
	ExecuteCommands("Tree_newick." + branchNames[k] + ".t:=" + branchLengths[k] + ";");
}

LikelihoodFunction  LikStateModel = (filt_data, Tree_newick);
Optimize (paramValues, LikStateModel);
fprintf (stdout, LikStateModel);
fprintf ("{{ name }}_hyout.txt", LikStateModel);
