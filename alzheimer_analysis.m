%% Big Data in Biomedical Sciences Final Report Code
% Buddy Pair 18: Anna Serbina () and Kristine McLaughlin (2837098)
% Last modified May 28, 2024

%% Week 2
% Load the GWAS GENE data into table
TGWAS = readtable('AD_magma.genes.txt', 'FileType', 'text');

% Only want columns GENE, CHR, ZSTAT, SYMBOL, P
TGWASclean = TGWAS;
TGWASclean.STOP=[];
TGWASclean.START=[];
TGWASclean.NSNPS=[];
TGWASclean.NPARAM = [];
TGWASclean.N = [];

% Make a GWAS plot
figure,hold on
pos=[1:size(TGWASclean.P)];
plot(pos,-log10(TGWASclean.P(:)),'.')
odd= find(mod(TGWAS.CHR,2)==1);
plot(pos(odd),-log10(TGWASclean.P(odd)),'.r')

for ii=1:22
tmp=find(TGWASclean.CHR==ii);
XT(ii)=round((tmp(end)-tmp(1))./2) + tmp(1);
XTlabels{ii}=(['CH',num2str(ii)]);
end

xticks(double(XT));
xticklabels(XTlabels);
y = -log10(10.^-8);
line([1,pos(end)],[y,y],'LineStyle','--');

hold off

% Distribution of p-values
figure, histogram(-log10(TGWASclean.P(:)));

% Get data for symbol APOE
APOE_loc = find(contains(TGWASclean.SYMBOL, 'APOE'));
APOE_zstat = TGWASclean.ZSTAT(APOE_loc)

% Identify highest z score
max_zstat = max(TGWASclean.ZSTAT)
max_zstat_ind = find(TGWASclean.ZSTAT == max_zstat);
symbol_maxZstat = TGWASclean.SYMBOL(max_zstat_ind)

% Find genes that reach significance threshold
threshold = 10.^-8;
significant_genes_ind = find(TGWASclean.P <= threshold);

SIGNIFICANT_AD_GENES = TGWASclean.SYMBOL(significant_genes_ind);
significant_genes_tbl = TGWASclean(significant_genes_ind, :)
% sort by ascending p value
significant_genes_tbl = sortrows(significant_genes_tbl, {'P'});

% Store results in text file
% writetable(significant_genes_tbl, "significant_genes.txt");

%% Week 3
% Load gene expression data into tables
TRAWGeneExpdata=readtable('ExpressionData.txt');
TMETAGeneExpdata=readtable('SamplesMetaData.txt');
Tprobes=readtable('ProbeNames.txt');

% Merging tables into one
% Put probid and corresponding gene_symbol into new table
tidx1 = find(contains(Tprobes.Properties.VariableNames,'probid'));
tidx2 = find(strcmp(Tprobes.Properties.VariableNames,'gene_symbol'));
tbl_probeNames = Tprobes(:,[tidx1, tidx2]);
clear tidx* Tprobes %do some clean up, we don't need these anymore

% Need column names of the key to be identical between tables
tidx1=find(contains(tbl_probeNames.Properties.VariableNames, "probid"));
tbl_probeNames.Properties.VariableNames{tidx1} = 'probe_id';
clear tidx*

% Merge probe names with expression data
tbl_expression_probeNames = join(tbl_probeNames, TRAWGeneExpdata, "Keys", "probe_id");

% Create table for one gene expression in particular brain region
RegionOfInterest='MTG';
meta_REGindx = (strcmp(TMETAGeneExpdata.regionDescriptions,RegionOfInterest));
tbl_meta_RegionOfInterest = TMETAGeneExpdata(meta_REGindx, :);

GeneOfInterest='APOE';
data_GENindx = strcmp(tbl_expression_probeNames.gene_symbol,GeneOfInterest);

% Get samples that match gene
tmpmatch_tissueSamples_indx = find(ismember(tbl_expression_probeNames.Properties.VariableNames, tbl_meta_RegionOfInterest.tissueSampleDescriptions));
tmpSampleVals = tbl_expression_probeNames(data_GENindx,tmpmatch_tissueSamples_indx);
tmpSampleVals = mean(tmpSampleVals{:,:});
tmpT = table(tmpSampleVals', 'VariableNames', {'APOEvals'});
% Add mean values into table
tbl_data_GEN_RegionOfInterest=[tbl_meta_RegionOfInterest, tmpT]

% Test for difference between controls and patients
% select CON and PAT values from our table
meta_CONindx = (strcmp(tbl_data_GEN_RegionOfInterest.groupID,'control'));
meta_PATindx = (strcmp(tbl_data_GEN_RegionOfInterest.groupID,'affected'));

% select CON and PAT values from our table
ConVal=tbl_data_GEN_RegionOfInterest.APOEvals(meta_CONindx,:);
PatVal=tbl_data_GEN_RegionOfInterest.APOEvals(meta_PATindx,:);

% APOE Histogram for each group
figure,
h1= histogram(ConVal);
hold on
h2= histogram(PatVal, 'BinEdges', h1.BinEdges); 

% Creating histogram for SORL1 gene in region PC
RegionOfInterest='PC';
meta_REGindx = (strcmp(TMETAGeneExpdata.regionDescriptions,RegionOfInterest));
tbl_meta_RegionOfInterest = TMETAGeneExpdata(meta_REGindx, :);

GeneOfInterest='SORL1';
data_GENindx = strcmp(tbl_expression_probeNames.gene_symbol,GeneOfInterest);

% Get samples that match gene
tmpmatch_tissueSamples_indx = find(ismember(tbl_expression_probeNames.Properties.VariableNames, tbl_meta_RegionOfInterest.tissueSampleDescriptions));
tmpSampleVals = tbl_expression_probeNames(data_GENindx,tmpmatch_tissueSamples_indx);
tmpSampleVals = mean(tmpSampleVals{:,:});
tmpT = table(tmpSampleVals', 'VariableNames', {'SORL1vals'});
% Add mean values into table
tbl_data_GEN_RegionOfInterest=[tbl_meta_RegionOfInterest, tmpT];

% Test for difference between controls and patients
% select CON and PAT values from our table
meta_CONindx = (strcmp(tbl_data_GEN_RegionOfInterest.groupID,'control'));
meta_PATindx = (strcmp(tbl_data_GEN_RegionOfInterest.groupID,'affected'));

% select CON and PAT values from our table
ConVal=tbl_data_GEN_RegionOfInterest.SORL1vals(meta_CONindx,:);
PatVal=tbl_data_GEN_RegionOfInterest.SORL1vals(meta_PATindx,:);

% SORL1 Histogram for each group
figure,
h1= histogram(ConVal, 'BinWidth', 500, 'FaceColor', "#1E88E5");
hold on
h2= histogram(PatVal, 'BinWidth', 500, 'FaceColor', "#D81B60");
legend('Control','Patient')
xlabel('SORL1 Expression Values')
ylabel('Frequency')
title('SORL1 Expression in Piriform Cortex between Controls and AD Patients')

% Barplot
x = 1:2;
data=[mean(ConVal), mean(PatVal)];
errlow=[mean(ConVal)-std(ConVal) mean(PatVal)-std(PatVal) ];
errhigh=[mean(ConVal)-std(ConVal) mean(PatVal)-std(PatVal) ];
figure, bar(x, data)
hold on
er=errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

% Perform a t-test between the two groups
[h,p,ci,stats]=ttest2(ConVal, PatVal)

%% Testing for all significant genes from GWAS

% Run on all significant genes in MTG
TABLE_GENE_EXPRESSION_SIGGENES_MTG = table();

for i = 1:length(SIGNIFICANT_AD_GENES)
    gene = SIGNIFICANT_AD_GENES{i};

    %zstat for a gene
    gene_loc = find(strcmp(TGWASclean.SYMBOL, gene));
    gene_zstat = TGWASclean.ZSTAT(gene_loc);
    
    %if there is no expression for a gene, skip
    if isempty(gene_zstat)
        continue
    end
    
    output_struct = select_gene_expression_values('MTG', gene, TMETAGeneExpdata, tbl_expression_probeNames);
    output_struct.Z_gwas = gene_zstat;
    output_tbl = struct2table(output_struct);

    TABLE_GENE_EXPRESSION_SIGGENES_MTG = [TABLE_GENE_EXPRESSION_SIGGENES_MTG; output_tbl];
end

% Run on all significant genes in every brain region
TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS = table();
all_regions = unique(TMETAGeneExpdata.regionDescriptions);

for region_idx = 1:length(all_regions)
    for gene_idx = 1:length(SIGNIFICANT_AD_GENES)
        region = all_regions{region_idx};
        gene = SIGNIFICANT_AD_GENES{gene_idx};

        % get zstat for a gene
        gene_loc = find(strcmp(TGWASclean.SYMBOL, gene));
        gene_zstat = TGWASclean.ZSTAT(gene_loc);
        
        % if there is no expression for a gene, skip
        if isempty(gene_zstat)
            continue
        end
        
        % get t-test value
        output_struct = select_gene_expression_values(region, gene, TMETAGeneExpdata, tbl_expression_probeNames);
        output_struct.Z_gwas = gene_zstat;
        output_tbl = struct2table(output_struct);
    
        % append to master table
        TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS = [TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS; output_tbl];
    end
end

num_tests = height(TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS);

% Filter out rows with missing values from the table
TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS = rmmissing(TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS);

% Filter out genes that show a significant level of expression, using p val
p_threshold = 0.05 / num_tests;
significant_ind = find(TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS.P_val <= p_threshold);

TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS_SIGNIF = TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS(significant_ind, :);
TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS_SIGNIF = sortrows(TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS_SIGNIF, {'P_val'})

% Get most significant gene for each region
TABLE_GENE_EXPRESSION_SIGGENES_SIGNI_largest = table();
for region_idx = 1:length(all_regions)
    % get rows of the desired region
    region_indices = find(strcmp(TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS.REGION, all_regions{region_idx}));
    region_tbl = TABLE_GENE_EXPRESSION_SIGGENES_ALLREGIONS(region_indices, :);
    
    % sort to find lowest p-value = most significant
    region_tbl = sortrows(region_tbl, {'P_val'});

    % add to TABLE_GENE_EXPRESSION_SIGGENES_SIGNI_largest
    TABLE_GENE_EXPRESSION_SIGGENES_SIGNI_largest = [TABLE_GENE_EXPRESSION_SIGGENES_SIGNI_largest; region_tbl(1,:)];
end

TABLE_GENE_EXPRESSION_SIGGENES_SIGNI_largest = sortrows(TABLE_GENE_EXPRESSION_SIGGENES_SIGNI_largest, {'P_val'});

%% Week 4
% Find top 10 most affecting genes in HIP region
RegionOfInterest='HIP';
meta_REGindx = (strcmp(TMETAGeneExpdata.regionDescriptions,RegionOfInterest));
tbl_meta_RegionOfInterest = TMETAGeneExpdata(meta_REGindx, :)

% select the samples out of the total set of probes and samples
tmpmatch_tissueSamples_indx = find(ismember(tbl_expression_probeNames.Properties.VariableNames, tbl_meta_RegionOfInterest.tissueSampleDescriptions));
tmpSampleVals = tbl_expression_probeNames(:, tmpmatch_tissueSamples_indx);

% define which samples are from the controls and which are from the patients
meta_CONindx = strcmp(tbl_meta_RegionOfInterest.groupID, 'control');
meta_PATindx = strcmp(tbl_meta_RegionOfInterest.groupID,'affected');

% select the values of these samples and perform a t-test between the two
tmpSampleVals_val=tmpSampleVals{:,:};
ConVal=tmpSampleVals_val(:,meta_CONindx);
PatVal=tmpSampleVals_val(:,meta_PATindx);
[a,b,c,d]=ttest2(ConVal', PatVal');

% sort the p-values of all the performed t-test, select the 10 lowest (so strongest effects)
[aa,bb] = sort(b);
tmptop10_probes = bb(1:10);

% List the genes that match the top10_probes
GENElist_top10 = tbl_probeNames(tmptop10_probes, 'gene_symbol');

%find the genes in the GWAS file and get their Zscore
clear tmpzstats
for ii=1:10
[a] = contains(TGWASclean.SYMBOL, GENElist_top10.gene_symbol{ii});
tmpzstats(ii)=mean(TGWASclean.ZSTAT(find(a)));
end

%average the values to get a mean ‘GWAS ZSCORE'; we take the absolute
nanmean(abs(tmpzstats))

%% Analyzing PPI
TPPI_intact = readtable('intact_int.txt');
TENGS_SYMBOL = readtable('ENGS_SYMBOL_list_intact.txt');
GENE_of_INTEREST = 'APOE';
clear PPIs_of_interest;
tel = 0;

% Find all indices where TENGS_SYMBOL.gene_name matches GENE_of_INTEREST
id = find(strcmp(TENGS_SYMBOL.gene_name, GENE_of_INTEREST));

% Initialize an empty array to store indices of matching ENGS
d = [];

% Iterate over each index in id
for i = 1:numel(id)
    % Find the indices where TPPI_intact.ensg1 matches the current ENGS
    id2 = find(strcmp(TPPI_intact.ensg1, TENGS_SYMBOL.ensg(id(i))));
    
    % Find the indices where TENGS_SYMBOL.ensg matches TPPI_intact.ensg2 using logical indexing
    d = [d; find(ismember(TENGS_SYMBOL.ensg, TPPI_intact.ensg2(id2)))];
end

for i = 1:numel(id)
    % Find the indices where TPPI_intact.ensg2 matches the current ENGS
    id3 = find(strcmp(TPPI_intact.ensg2, TENGS_SYMBOL.ensg(id(i))));
    
    % Find the indices where TENGS_SYMBOL.ensg matches TPPI_intact.ensg2 using logical indexing
    d = [d; find(ismember(TENGS_SYMBOL.ensg, TPPI_intact.ensg1(id3)))];
end

% Display the corresponding gene names
TENGS_SYMBOL.gene_name(d)

%% Find PPIs for all genes from Week 1 GWAS
clear PPIs_of_interest; tel=0;
for jj = 1:length(SIGNIFICANT_AD_GENES)
    PROT_of_INTEREST = SIGNIFICANT_AD_GENES(jj);
    id=find(strcmp(TENGS_SYMBOL.gene_name, PROT_of_INTEREST));
    if nnz(id)==0; continue; end

    id1=find(ismember(TPPI_intact.ensg1, TENGS_SYMBOL.ensg(id(i))));
    d1=find(ismember(TENGS_SYMBOL.ensg, TPPI_intact.ensg2(id1)));

    id2=find(ismember(TPPI_intact.ensg2, TENGS_SYMBOL.ensg(id(i))));
    d2=find(ismember(TENGS_SYMBOL.ensg, TPPI_intact.ensg1(id2)));

    d=[d1;d2];
    d=unique(d);

    for ii=1:size(d,1)
        tel=tel+1;
        PPIs_of_interest(tel) = TENGS_SYMBOL.gene_name(d(ii));
    end
end
PPIs_of_interest = unique(PPIs_of_interest);

% put the data into a table
TPPI = table(PPIs_of_interest(:));
TPPI.Properties.VariableNames = {'gene_name'};

% Test whether the proteins show alterned expression levels
RegionOfInterest='MTG';
for ii=1:size(PPIs_of_interest,2)
    a=find(strcmp(tbl_expression_probeNames.gene_symbol,PPIs_of_interest{ii}));
    if isempty(a);continue;end
    meta_REGindx = (strcmp(TMETAGeneExpdata.regionDescriptions, RegionOfInterest));
    tbl_meta_RegionOfInterest = TMETAGeneExpdata(meta_REGindx,:);
    tmpmatch_tissueSamples_indx = find(ismember(tbl_expression_probeNames.Properties.VariableNames, tbl_meta_RegionOfInterest.tissueSampleDescriptions));
    tmpSampleVals = tbl_expression_probeNames(a, tmpmatch_tissueSamples_indx);
    meta_CONindx = (strcmp(tbl_meta_RegionOfInterest.groupID,'control'));
    meta_PATindx = (strcmp(tbl_meta_RegionOfInterest.groupID, 'affected'));
    ConVal=mean(tmpSampleVals{:,meta_CONindx},1);
    PatVal=mean(tmpSampleVals{:,meta_PATindx},1);
    [h,p,ci,stats]=ttest2(ConVal, PatVal);
    TPPI.GENexpdstats(ii) = stats.tstat;
    clear tmp*
end

TPPI = sortrows(TPPI, {'GENexpdstats'});

%% 4.3 PPI Networks
% get unique proteins particularly thought to be involved in AD
TPPI_AD_intact = readtable('alz_intact_int.txt');
ADproteins = unique([TPPI_AD_intact.ensg1, TPPI_AD_intact.ensg2]);

% build network
unique_ensg1 = unique(TPPI_AD_intact.ensg1);
unique_ensg2 = unique(TPPI_AD_intact.ensg2);
unique_ensg_all = unique([unique_ensg1; unique_ensg2]);

PPI_AD_matrix = zeros(size(unique_ensg_all,1));

% make matrix bidirectional
for ii=1:size(TPPI_AD_intact,1)
    t1=TPPI_AD_intact.ensg1{ii};
    t2=TPPI_AD_intact.ensg2{ii};
    PPI_AD_matrix(find(contains(unique_ensg_all,t1)), find(contains(unique_ensg_all,t2)))=TPPI_AD_intact.score(ii);
end
PPI_AD_matrix_bidirectional = triu(PPI_AD_matrix, 0) + triu(PPI_AD_matrix)';

% make matrix binary - whether there is a PPI or not
PPI_AD_matrix_bidirectional_binary = double(PPI_AD_matrix_bidirectional>0);

%%
G = graph(PPI_AD_matrix_bidirectional_binary);
% G = graph(PPI_AD_matrix_bidirectional);
G.Nodes.Names = unique_ensg_all;
figure;
plot(G,'NodeLabel', G.Nodes.Names, 'LineWidth', G.Edges.Weight,'Layout', 'force')

% Calculate the degree of each node in the graph
node_degrees = degree(G);

% Find the node with the highest degree
[max_degree, max_degree_node] = max(node_degrees);

% Display the result
fprintf('Node with the most edges: %s (Degree: %d)\n', G.Nodes.Names{max_degree_node}, max_degree);
%ENSG00000142192 gene

%% Explore weights of PPIs
allPPIvalues = PPI_AD_matrix(PPI_AD_matrix ~= 0);

% Find the maximum PPI value in the matrix
[max_PPI_value, max_PPI_index] = max(PPI_AD_matrix_bidirectional(:));

% Convert linear index to row and column indices
[max_row, max_col] = ind2sub(size(PPI_AD_matrix_bidirectional), max_PPI_index);

% Retrieve the corresponding proteins using the row and column indices from the unique_ensg_all array
protein1 = unique_ensg_all(max_row)
protein2 = unique_ensg_all(max_col)

%ENSG00000142192 has the strongest connections

figure, histogram(allPPIvalues)

[a,b]=modularity_und(PPI_AD_matrix_bidirectional);

% a is optimal community structure, b is maximized modularity
num_modules = max(a)

% Create a cell array with protein names and their corresponding Ci values
protein_Ci_pairs = [G.Nodes.Names(:), num2cell(a(:))];

% Sort the cell array based on the Ci values
nodes_grouped_together = sortrows(protein_Ci_pairs, 2);

%% Week 5: Integration with Brain Gene Expression
% Load transcriptomics data
AHBA = load('AHBA_transcriptomics_114atlas.mat');

% Only interested in left hemisphere data bc complete
lh_indices = find(contains(AHBA.regionDescription, "lh"));
AHBA_lefthemi = AHBA.gene_expression_region_FS_z(lh_indices, :, :);
AHBA_lefthemi_group_average = nanmean(AHBA_lefthemi,3);

% Make histogram with APOE gene
APOEGENE = find(strcmp(AHBA.gene_symbol, 'APOE'));
APOE_expression_DK114_lh = AHBA_lefthemi_group_average(:,APOEGENE);

figure, histogram(APOE_expression_DK114_lh)

colormap = [ones(1001,1), linspace(0.7529, 0.2510, 1001)', ones(1001,1)]; %these are RGB colors - you can play around with this making your figure look nicer and nicer!
regionNames = AHBA.regionDescription; %list of names of the regions
corticalMap_LR = [APOE_expression_DK114_lh; zeros(57,1)]'; %the script is for both hemispheres and thus expects 114 region values.
plotSurfaceBoth(regionNames, corticalMap_LR, colormap)

%I think the brain images show abnormal transcription levels in cerebellum
%which agrees with the info on the Internet

% Extract top 10 regions where APOE is most expressed in cortex
[aa,bb] = sort(APOE_expression_DK114_lh,'descend');
expressionLevels = aa(1:10);
APOE_high_expr_ind = bb(1:10); 
ABHA_lh_regionDescr = AHBA.regionDescription(lh_indices);
APOE_high_expr_regions = ABHA_lh_regionDescr(APOE_high_expr_ind);
rank = (1:10)';
APOE_high_expr_regions_table = table(rank, APOE_high_expr_regions, expressionLevels, 'VariableNames', {'Rank', 'RegionName', 'ExpressionLevel'});

%% Generate high region expression levels for other genes from GWAS
%5.1.11
% Initialize variables to store expression levels and final table
SIGGENES_topABHAexprregion = table();

for jj = 1:length(SIGNIFICANT_AD_GENES)
    PROT_of_INTEREST = SIGNIFICANT_AD_GENES{jj}; % Note the curly braces for cell array indexing
    GENE = find(strcmp(AHBA.gene_symbol, PROT_of_INTEREST));
    
    if ~isempty(GENE) % Check if GENE is not empty
        GENE_expression_DK114_lh = AHBA_lefthemi_group_average(:, GENE);
        [expressionLevel, ind] = sort(GENE_expression_DK114_lh, 'descend');
        ind = ind(1);
        topRegion = AHBA.regionDescription(ind);
        
        % Create table for current gene and append to final_table
        geneTable = table({PROT_of_INTEREST}, topRegion, 'VariableNames', {'GENESYMBOL', 'TOPREGION_AHBAEXPRESSION'});
        SIGGENES_topABHAexprregion = [SIGGENES_topABHAexprregion; geneTable];
    end
end

% Display or use the final table as needed
disp(SIGGENES_topABHAexprregion);

%% 5.2 VISUALISATION
% Create table for one gene expression in particular brain region
RegionOfInterest='MTG';
meta_REGindx = (strcmp(TMETAGeneExpdata.regionDescriptions,RegionOfInterest));
tbl_meta_RegionOfInterest = TMETAGeneExpdata(meta_REGindx, :);

GeneOfInterest='APOE';
data_GENindx = strcmp(tbl_expression_probeNames.gene_symbol,GeneOfInterest);

% Get samples that match gene
tmpmatch_tissueSamples_indx = find(ismember(tbl_expression_probeNames.Properties.VariableNames, tbl_meta_RegionOfInterest.tissueSampleDescriptions));
tmpSampleVals = tbl_expression_probeNames(data_GENindx,tmpmatch_tissueSamples_indx);
tmpSampleVals = mean(tmpSampleVals{:,:});
tmpT = table(tmpSampleVals', 'VariableNames', {'APOEvals'});
% Add mean values into table
tbl_data_GEN_RegionOfInterest=[tbl_meta_RegionOfInterest, tmpT];

% Test for difference between controls and patients
% select CON and PAT values from our table
meta_CONindx = (strcmp(tbl_data_GEN_RegionOfInterest.groupID,'control'));
meta_PATindx = (strcmp(tbl_data_GEN_RegionOfInterest.groupID,'affected'));

%select CON and PAT values from our table
ConVal=tbl_data_GEN_RegionOfInterest.APOEvals(meta_CONindx,:);
PatVal=tbl_data_GEN_RegionOfInterest.APOEvals(meta_PATindx,:);
plotdata.Controls = ConVal;
plotdata.Patients = PatVal;
figure; violinplot(plotdata);
ylabel('Average expression level');
title(sprintf('%s gene expression in the %s', GeneOfInterest, RegionOfInterest));


%%
% 5.3 Examining if significant changes in gene expression in 3 cell types
% Load in data files, rename column names bc need unique labels
AST_DATA = readtable('GSE102956_COUNTS_AST.txt', 'ReadVariableNames', false);
AST_DATA.Properties.VariableNames = {'GENE', 'AST_E3_1', 'AST_E3_2', 'AST_E3_3', 'AST_E4_1', 'AST_E4_2', 'AST_E4_3'};
 
NEU_DATA = readtable('GSE102956_COUNTS_NEU.txt', 'ReadVariableNames', false);
NEU_DATA.Properties.VariableNames = {'GENE', 'NEU_E3_1', 'NEU_E3_2', 'NEU_E3_3', 'NEU_E4_1', 'NEU_E4_2', 'NEU_E4_3'};

MIC_DATA = readtable('GSE102956_COUNTS_MIC.txt', 'ReadVariableNames', false);
MIC_DATA.Properties.VariableNames = {'GENE', 'MIC_E3_1', 'MIC_E3_2', 'MIC_E3_3', 'MIC_E3_4', 'MIC_E4_1', 'MIC_E4_2', 'MIC_E4_3', 'MIC_E4_4'};

% Generate table with significant different gene expressions for each cell
% type
iPSCsummaryTbl = table();
for cell_type = ["AST", "NEU", "MIC"]

    if strcmp(cell_type, "AST")
        iPSC_DATA = AST_DATA;
    elseif strcmp(cell_type, "NEU")
        iPSC_DATA = NEU_DATA;
    else
        iPSC_DATA = MIC_DATA;
    end

    % Extract columns for APOE3 and APOE4
    APOE3neurons = contains(iPSC_DATA.Properties.VariableNames, 'E3');
    APOE3columnsInTable = iPSC_DATA.Properties.VariableNames(APOE3neurons);
    APOE4neurons = contains(iPSC_DATA.Properties.VariableNames, 'E4');
    APOE4columnsInTable = iPSC_DATA.Properties.VariableNames(APOE4neurons);
    
    Data3=iPSC_DATA{:,APOE3columnsInTable};
    Data4=iPSC_DATA{:,APOE4columnsInTable};
    
    % Remove rows with any zeros
    rowsWith03 = any(Data3(:) == 0, 2);
    rowsWith04 = any(Data4(:) == 0, 2);
    
    Data3 = Data3(~rowsWith03);
    Data4 = Data4(~rowsWith04);
    
    % Perform t test between APOE3 and APOE4 cell lines
    [H, P, CI, STATS] = ttest2(Data3', Data4');
    
    % Get significant genes
    alpha = 0.05/nnz(P); % Bonferroni correction
    iPSCdiffExpGenes = (iPSC_DATA.GENE(P<alpha));
    iPSCdiffExpGenes_SIGN = sign(STATS.tstat(P<alpha));

    if isempty(iPSCdiffExpGenes)
        continue;
    end
    
    % Add other statistics to a table
    oneCellSummaryTbl = cell2table(iPSCdiffExpGenes);
    oneCellSummaryTbl.Properties.VariableNames(1) = "Gene"; % rename column
    oneCellSummaryTbl.P_value = P(P<alpha).'; % add P values
    oneCellSummaryTbl.Cell_type = repmat(cell_type, size(oneCellSummaryTbl,1), 1);

    % Add this table to complete summary table
    iPSCsummaryTbl = [iPSCsummaryTbl; oneCellSummaryTbl];
end

disp(iPSCsummaryTbl);

%% WEEK 6
% Load MRI data
MRIdata=load('SURASIS_grouped_regionProperties_FS_aparc2.mat'); % raw data
tbl_meta_data_mri=readtable('SURASIS_mr_sessions.csv'); % MRI sessions of subjects
tbl_meta_data_clinical=readtable('SURASIS_clinical_data.csv'); % clinical data of subjects

% Only subjects of which we have MRI data for
[a,b]=ismember(MRIdata.subjects, tbl_meta_data_mri.MRID);
tbl_meta_data_mri_filtered = tbl_meta_data_mri(b, :);

% find for all subjects in MRIdata their PATIENT CONTROL STATUS and add this
% to the meta_data_mri_sessions table
for s = 1: numel(MRIdata.subjects)
    tmpsub_mri_session = MRIdata.subjects{s};
    % step 1: find the index of tmpsub_mri_session in tbl_meta_data_mri_filtered. Call it for example tmp1
    tmp1=find(strcmp(tbl_meta_data_mri_filtered.MRID, tmpsub_mri_session));
    
    % step 2: find the subcode (call it for example tmpsubcode) that matches tmp1
    tmpsubcode = tbl_meta_data_mri_filtered.Subject{tmp1};
    
    % step 3: find the clinical status of this subject(tmpsubcode) in tbl_meta_data_clinical (where is this information in this table, i.e. which column to use?)
    tmp2=find(strcmp(tbl_meta_data_clinical.Subject, tmpsubcode));
    
    % step 4: extract the status, call it for example tmpstatus
    tmpstatus=tbl_meta_data_clinical.dx1(tmp2);
    
    % We can simply take the first of the measured status to determine
    % whether one is a control or a patient
    tmpstatus=tmpstatus{1};
    
    if strcmp(tmpstatus, 'Cognitively normal')
        tbl_meta_data_mri_filtered.status{tmp1} ='CON';
    elseif strcmp(tmpstatus, 'AD Dementia')
        tbl_meta_data_mri_filtered.status{tmp1} = 'AD';
    else
        tbl_meta_data_mri_filtered.status{tmp1} = 'Other';
    end
end

%% Exploring if MTG region shows structural differences in patients
% find the region middle temporal gyrus
regions = find(contains(MRIdata.regionDescriptions,'middletemporal'));
display(regions)
dataProp = find(contains(MRIdata.propertyDescriptions,'GrayVol'));
% there is one region for the left and one for the right hemisphere
ROIvol = squeeze(MRIdata.regionProperties(regions,dataProp,:));
ROIvol = mean(ROIvol);

% remove outliers 3 standard deviations from the mean
outliers = (ROIvol < nanmean(ROIvol) - 3*nanstd(ROIvol) | ROIvol > nanmean(ROIvol) + 3*nanstd(ROIvol));
ROIvol(outliers) = NaN; % we exclude these values from further analysis

% determine who are Controls and who are the Patients
CONindx = find(strcmp(tbl_meta_data_mri_filtered.status,'CON'));
PATindx = find(strcmp(tbl_meta_data_mri_filtered.status,'AD'));

plotdata.Controls = ROIvol(CONindx);
plotdata.Patients = ROIvol(PATindx);
figure; violinplot(plotdata);
ylabel('Volume (mm3)');
title('MTG Regional gray matter volume');

% difference in volume
[h,p,~,stats]=ttest2(ROIvol(CONindx), ROIvol(PATindx));

%%
% difference in volume for all regions

% Extract the regional properties for all regions and all subjects
ROIvol = squeeze(MRIdata.regionProperties(:, dataProp, :));

% Detect outliers
outliers = (ROIvol < nanmean(ROIvol, 2) - 3*nanstd(ROIvol, 0, 2)) | (ROIvol > nanmean(ROIvol, 2) + 3*nanstd(ROIvol, 0, 2));

% Remove outliers from the control and patient indices
CONindx(find(ismember(CONindx, find(outliers)))) = [];
PATindx(find(ismember(PATindx, find(outliers)))) = [];

for regionIdx = 1:length(MRIdata.regionDescriptions)
    [~, ~, ~, stats] = ttest2(ROIvol(regionIdx, CONindx), ROIvol(regionIdx, PATindx));
    tStats(regionIdx) = stats.tstat;
end

% Create a table with region names and t-test statistics
AD_pathology_score_table = table(MRIdata.regionDescriptions, tStats', 'VariableNames', {'Region', 'TTestScore'});

% Get top 10 regions with largest volume differences
[aa,bb] = sort(AD_pathology_score_table.TTestScore,'descend');
top_10_regions = AD_pathology_score_table(bb(1:10), :)

cm=[0:0.001:1; 0:0.001:1 ;0:0.001:1]'; cm=flipud(cm); % black and white
%cm = [ones(1001,1), linspace(0.7529, 0.2510, 1001)', ones(1001,1)]; % pink
plotSurfaceBoth_DK68(MRIdata.regionDescriptions, tStats, cm)


%% Exploring second q: Is the expression pattern of APOE across cortical 
% regions associated with the pattern of AD structural changes across the cortex?
figure,histogram(APOE_expression_DK114_lh)

% map the 68 regions to the 114 regions we have AHB expr data for
APOE_expression_DK68_lh = [];
for ii=1:numel(MRIdata.regionDescriptions)/2 %changed size to numel
    DK68region = MRIdata.regionDescriptions{ii};
    AHBAregions = find(contains(AHBA.regionDescription, DK68region));
    APOE_expression_DK68_lh(ii) = nanmean(APOE_expression_DK114_lh(AHBAregions));
end

AD_pathology_score_DK68 = tStats;
% Average atrophy pattern across the two hemispheres left and right to get 1 vector
lhregions_idx=find(contains(MRIdata.regionDescriptions, 'lh-'));
rhregions_idx=find(contains(MRIdata.regionDescriptions, 'rh-'));
AD_pathology_score_DK68_lhrh = nanmean([AD_pathology_score_DK68(lhregions_idx);AD_pathology_score_DK68(rhregions_idx)]);
figure, plot(APOE_expression_DK68_lh',AD_pathology_score_DK68_lhrh,'.');

%% WEEK 6 PART 2 Visualization
f = figure;
hold on

x=APOE_expression_DK68_lh';
y=AD_pathology_score_DK68_lhrh';
s = scatter(x,y);
l = lsline;
f.Color = 'white';

color1 = [72 150 236] ./ 255;
s.MarkerFaceColor = color1;
s.MarkerEdgeColor = color1;
s.MarkerEdgeAlpha = .3;
s.MarkerFaceAlpha = .3;
l.Color = 'black';
l.LineWidth = 1;

% Add labels to the x-axis
xlabel('APOE expression');
ylabel('AD pathology score');

a = gca;
a.XAxis.LineWidth = 1;
a.YAxis.LineWidth = 1;
a.FontSize = 12;

% We need to sort the values of x in ascending order for the
%function patch() below - this has no effect on our fit function.
[~, I] = sort(x);
xsorted = x(I);
ysorted = y(I);
fitresult = fit(xsorted, ysorted, 'poly1');
p95 = predint(fitresult,xsorted,0.95,'observation','off');
p90 = predint(fitresult,xsorted,0.90,'observation','off');

plot(xsorted,p95,'r--')
plot(xsorted,p90,'r--')

patch([xsorted; flipud(xsorted)], [p90(:,1); flipud(p90(:,2))], 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.05);
patch([xsorted; flipud(xsorted)], [p95(:,1); flipud(p95(:,2))], 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.1);
% replot the points to bring them on top
plot(xsorted,ysorted,'.b')
title('Correlation Between APOE Expression and Region Vulnerability in AD Pathology');


%% Part 3 Statistics
s = regstats(APOE_expression_DK68_lh, AD_pathology_score_DK68_lhrh);
disp(s.tstat.pval)

%% Part 4 CONNECTOMICS
% there is no real quick way to avoid hardcoding here (except if you
%have Linux or Mac), so lets ‘hard’ define the number of APOE3 and
%APOE4 carriers. These values will never change so this is justified to do.
nAPOE3_carriers=30;
nAPOE4_carriers=25;

% load in the META data on AGE
tmp1 = readtable('ages_list.txt','ReadVariableNames',false);
tmp1.Properties.VariableNames{1}='Age';
% load in the META data on GENDER
tmp2 = readtable('genders_list.txt','ReadVariableNames',false);
tmp2.Properties.VariableNames{1}='Gender';

% make a binary variable for GENDER, FEMALE=1, MALE=0
tmp3=tmp2;
gen=double(strcmp(tmp2.Gender,'Female'));
tmp3.Gender=gen;
tmp3.Properties.VariableNames{1}='GenderBin';
%put all demographics in one table
tblDemoGraphics = [tmp1,tmp2, tmp3];
clear tmp*

%split the Demographics table into 2 tables 1 for APOE3 and 1 for
% APOE4. We can also make 1 table with patient status in it as an
% alternative
tblDemoGraphics_APOE3 = tblDemoGraphics(1:nAPOE3_carriers,:);
tblDemoGraphics_APOE4 = tblDemoGraphics(nAPOE3_carriers+1:end,:);
%load in the META data on the RegionNames
regionNames = readtable('APOE-4_9_region_names_full_file.txt','ReadVariableNames', false);
regionNames_cleaned=[];
for ii=1:size(regionNames,1)
    tmp=regionNames{ii,:};
    totaltmp=[];
    for jj=1:size(tmp,2)
        if isempty(tmp{jj});continue;end
        totaltmp=[totaltmp,tmp{jj},'_'];
    end
    %remove the last
    totaltmp(end)=[];
    %our regionNames
    regionNames_cleaned{ii}=totaltmp;
end

%put regionNames in a table and clean all intermediate steps
regionNames=regionNames_cleaned';
tbl_regionNames=cell2table(regionNames);
clear regionNames_cleaned *tmp*
%load in the connectivity matrices for the APOE3 and APOE4 group
for ii=1:nAPOE3_carriers
    tmpconnectivity_matrix=load(['APOE-3_',num2str(ii),'_connectivity_matrix_file.txt']);
    APOE3_matrices(ii,:,:)=tmpconnectivity_matrix;
end
for ii=1:nAPOE4_carriers
    tmpconnectivity_matrix=load(['APOE-4_',num2str(ii),'_connectivity_matrix_file.txt']); %hint almost the same as in the APOE3 carrier loop
    APOE4_matrices(ii,:,:)=tmpconnectivity_matrix;
end

% Examine difference in connectivity APOE carriers in hippocampus
RegionOfInterest='hippocampus';
idx_roi = (strfind(tbl_regionNames.regionNames, RegionOfInterest));
idx_roi = find(~cellfun(@isempty, idx_roi));

%get the connectivity of all APOE3_carriers of the left hippocampus
CONN_ROI_APOE4 = nansum(squeeze(nansum(squeeze(APOE4_matrices(:,idx_roi,:)),2)),2);
CONN_ROI_APOE3 = nansum(squeeze(nansum(squeeze(APOE3_matrices(:,idx_roi,:)),2)),2);

% regress out effect of gender
s=regstats(CONN_ROI_APOE4, tblDemoGraphics_APOE4.GenderBin);
CONN_ROI_APOE4 = s.r;
s=regstats(CONN_ROI_APOE3, tblDemoGraphics_APOE3.GenderBin);
CONN_ROI_APOE3 = s.r;

% plot results
figure,plot(tblDemoGraphics_APOE3.Age,CONN_ROI_APOE3,'.'); lsline
hold on
plot(tblDemoGraphics_APOE4.Age,CONN_ROI_APOE4,'.r'); lsline
title("Hippocampus Connectivity in APOE carriers/non-carriers");

% Repeat analysis for precuneous region
RegionOfInterest='precuneous';
idx_roi = (strfind(tbl_regionNames.regionNames, RegionOfInterest));
idx_roi = find(~cellfun(@isempty, idx_roi));

CONN_ROI_APOE4 = nansum(squeeze(nansum(squeeze(APOE4_matrices(:,idx_roi,:)),2)),2);
CONN_ROI_APOE3 = nansum(squeeze(nansum(squeeze(APOE3_matrices(:,idx_roi,:)),2)),2);

s=regstats(CONN_ROI_APOE4, tblDemoGraphics_APOE4.GenderBin);
CONN_ROI_APOE4 = s.r;
s=regstats(CONN_ROI_APOE3, tblDemoGraphics_APOE3.GenderBin);
CONN_ROI_APOE3 = s.r;

figure,plot(tblDemoGraphics_APOE3.Age,CONN_ROI_APOE3,'.'); lsline
hold on
plot(tblDemoGraphics_APOE4.Age,CONN_ROI_APOE4,'.r'); lsline
title("Precuneous Connectivity in APOE carriers/non-carriers");

%% WEEK 7
MRIdata = load('SURASIS_grouped_regionProperties_FS_aparc2.mat');
% Unsupervised learning - using k means clustering
Region1 = 'ctx-lh-medialorbitofrontal';
Region2 = 'ctx-lh-middletemporal';
idxRegion1=find(strcmp(MRIdata.regionDescriptions, Region1));
idxRegion2=find(strcmp(MRIdata.regionDescriptions, Region2));
Reg1val=squeeze(MRIdata.regionProperties(idxRegion1, dataProp, :));
Reg2val=squeeze(MRIdata.regionProperties(idxRegion2, dataProp, :));
Reg1val= (Reg1val - mean(Reg1val)) ./ std(Reg1val);
Reg2val= (Reg2val - mean(Reg2val)) ./ std(Reg2val);

% Remove outliers that are more than 3 stdev away from the mean
Reg1val_ind = abs(Reg1val) <= 3; 
Reg2val_ind = abs(Reg2val) <= 3;

% Find common valid indices
validIndices = Reg1val_ind & Reg2val_ind;

% Filter the values using the common valid indices
Reg1val_filtered = Reg1val(validIndices);
Reg2val_filtered = Reg2val(validIndices);

% K means clustering, generate plot
figure,plot(Reg1val_filtered, Reg2val_filtered,'.')
X=[Reg1val, Reg2val];
[idxcluster]=kmeans(X, 2);
figure,plot(Reg1val, Reg2val,'.')
hold on
plot(Reg1val((idxcluster==2)), Reg2val((idxcluster==2)),'.r')
plot(Reg1val((idxcluster==1)), Reg2val((idxcluster==1)),'.b')

%% SUPERVISED LEARNING

testsetsize=40;
trainsetsize=60;
tmp=randperm(length(CONindx));
trainsetCON=CONindx(tmp(1:trainsetsize));
testsetCON=CONindx(tmp(trainsetsize+1:testsetsize+trainsetsize));
tmp=randperm(length(PATindx));
trainsetPAT=PATindx(tmp(1:trainsetsize));
testsetPAT=PATindx(tmp(trainsetsize+1:testsetsize+trainsetsize));

ROIvol = squeeze(MRIdata.regionProperties(:,dataProp,:));
ROIvol = rescale(ROIvol);

% make an output vector Y telling who are patients and who are controls
Ytrain = zeros(length([trainsetPAT; trainsetCON]),1);
Ytrain(1:length(trainsetPAT))=1;
% make it a boolean (the fitclinear function wants it that way)
Ytrain=(Ytrain==1);
% set the input data X
Xtrain=[ROIvol(:,trainsetPAT), ROIvol(:,trainsetCON)];
Xtrain=Xtrain';
% do the same for the test set, call this Ytest and Xtest
Ytest = zeros(length([testsetPAT; testsetCON]),1);
Ytest(1:length(testsetPAT))=1;
Ytest=(Ytest==1);
Xtest=[ROIvol(:,testsetPAT), ROIvol(:,testsetCON)];
Xtest=Xtest';

[MDL, fitinfo] = fitclinear(Xtrain, Ytrain,'Regularization','ridge')
%%
MLD_train = predict(MDL,Xtrain);
%inssampleerror = 1-nnz((predict(MDL,Xtrain) == Ytrain))./length(Xtrain)
inssampleerror = loss(MDL, Xtrain , Ytrain)
% Number of correctly classified points in the test set
num_correct_train = nnz(MLD_train == Ytrain)

% Number of incorrectly classified points in the test set
num_incorrect_train = nnz(MLD_train ~= Ytrain)


outofsampleerror = loss(MDL, Xtest , Ytest) %lower than on the train set which is surprising
MLD_test = predict(MDL,Xtest);
% Number of correctly classified points in the test set
num_correct_test = nnz(MLD_test == Ytest)

% Number of incorrectly classified points in the test set
num_incorrect_test = nnz(MLD_test ~= Ytest)

%%
%NEURAL NETWORKS
% divide data into development/test/validation sets
trainRatio = 0.60;
valRatio = 0.20;
testRatio = 0.20;

% Create train/test/val set.
% Note we set the default seed to make sure we get the same train/test/val subjects every run.
fprintf('random seed set to default.\n');
rng('default');
% divide the subjects of the dataset into training/validation/testsets
nsubjects = height(tbl_meta_data_mri_filtered);
subjects_study.train = false(nsubjects, 1);
subjects_study.val = false(nsubjects, 1);
subjects_study.test = false(nsubjects, 1);
[xtrain, xval, xtest] = dividerand(1:nsubjects, trainRatio, ... 
    valRatio, testRatio);
subjects_study.train(xtrain) = true;
subjects_study.val(xval) = true;
subjects_study.test(xtest) = true;

dataProp = find(contains(MRIdata.propertyDescriptions,'GrayVol'));
ROIvol = MRIdata.regionProperties(:,dataProp,:);
% ROIvol = Array ( regions X properties X subjects )
train_data = ROIvol;
% We normalize the data
train_data = (train_data - mean(train_data, 3)) ./std(train_data, [], 3);
%%
train_data
n_regions = 68; %size(train_data, 1);
n_properties = 1; %size(train_data, 2);

% define network architecture:
solvername = 'sgdm';
layers = [imageInputLayer([n_regions n_properties 1], 'name', 'inputlayer')
fullyConnectedLayer(10, 'name', 'first layer')
reluLayer('name', 'relu1');
fullyConnectedLayer(2, 'name', 'second layer');
softmaxLayer('name', 'softmax');
classificationLayer('name', 'classification')];
% Set the options for the training
options = {'Verbose', true, ...
    'Plots','training-progress', ...
'shuffle', 'every-epoch'};

options_study = trainingOptions(solvername,...
options{:}, ...
'ValidationData', ...
{permute(train_data(:,:,subjects_study.val), [1 2 4 3]), ...
categorical(strcmp(tbl_meta_data_mri_filtered.status(subjects_study.val), 'AD'))});

[trained_net, training_results] = trainNetwork( ...
permute(train_data(:,:,subjects_study.train), [1 2 4 3]), ...
categorical(strcmp(tbl_meta_data_mri_filtered.status(subjects_study.train), 'AD')),...
layers, options_study);

%% Evaluate Model
XTest = permute(train_data(:,:,subjects_study.test), [1 2 4 3]);
YTest = categorical(strcmp(tbl_meta_data_mri_filtered.status(subjects_study.test), 'AD'));
[YPredicted,probs] = classify(trained_net, XTest);
% testError = 1 - accuracy
testError = 1 - mean(YPredicted == YTest)
accuracy = 1 - testError
confusionmat(YTest, YPredicted)

%%
% use the function make_fair which balances the number of controls and patients during the training phase.
% Ensure same number of patients and controls in the val and train set (NOT necessary in test set).
fair_val = make_fair(strcmp(tbl_meta_data_mri_filtered.status,'AD'), subjects_study.val);
fair_train = make_fair(strcmp(tbl_meta_data_mri_filtered.status, 'AD'), subjects_study.train);
% Make 3 lists of CON and CASES (i.e. Patients) out of this
CON_CASE_list_fairval = strcmp(tbl_meta_data_mri_filtered.status(fair_val), 'AD');
CON_CASE_list_fairtrain = strcmp(tbl_meta_data_mri_filtered.status(fair_train), 'AD');
CON_CASE_list_fairtest = strcmp(tbl_meta_data_mri_filtered.status(subjects_study.test), 'AD');

% You can try out some parameter settings by including or excludingthem:
% (help trainingOptions)
options = {'Verbose', true, ...
'Plots','training-progress', ...
'shuffle', 'every-epoch', ...
'InitialLearnRate', 0.1, ...
'L2Regularization', 0.03, ...
'Momentum', 0.95, ...
'ValidationPatience', 30};

options_study = trainingOptions(solvername,...
options{:}, ...
'ValidationData', ...
{permute(train_data(:,:,fair_val), [1 2 4 3]), ...
categorical(CON_CASE_list_fairval)});

[trained_net, training_results] = trainNetwork( ...
permute(train_data(:,:,fair_train), [1 2 4 3]), ...
categorical(CON_CASE_list_fairtrain), ...
layers, options_study);
% Evaluate again!
XTest = permute(train_data(:,:,subjects_study.test), [1 2 4 3]);
YTest = categorical(CON_CASE_list_fairtest);
[YPredicted,probs] = classify(trained_net, XTest);
testError = 1 - mean(YPredicted == YTest);
accuracy = 1 - testError
confusionmat(YTest, YPredicted)
