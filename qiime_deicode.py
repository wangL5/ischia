#5/22/20
#Using DEICODE in Qiime2 to calculate Robust Aitchison beta diversity

#Required files:
#Sample metadata: shoot_samples_keep.txt
#Count_table.qza: id_shoot_countfiltered_table.qza
#Taxonomy: taxonomy.qza

#1. Convert DADA2/phyloseq objects to qiime2 objects:
#https://forum.qiime2.org/t/importing-dada2-and-phyloseq-objects-to-qiime-2/4683
 
#2. Download DEICODE qiime2 plugin
#https://library.qiime2.org/plugins/deicode/19/
 
#3. Within qiime2 environment:
qiime deicode rpca \
    --i-table shoot-countfiltered-table.qza \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot ischia_shoot_ordination.qza \
    --o-distance-matrix ischia_shoot_distance.qza
 
#output
#Saved PCoAResults % Properties('biplot') to: ischia_shoot_ordination.qza
#Saved DistanceMatrix to: ischia_shoot_distance.qza
 
#biplot visualization
 
qiime emperor biplot \
    --i-biplot ischia_shoot_ordination.qza \
    --m-sample-metadata-file shoot_samples_metadata.txt \
    --m-feature-metadata-file taxonomy.qza \
    --o-visualization ischia_biplot.qzv \
    --p-number-of-features 5
 
#number of features: number of most important features (arrows) to display in the ordination. "Importance" is calculated for each feature based on the vector's magnitude, i.e. euclidean distance from origin. default 5.
 
#output
#Saved Visualization to: ischia_biplot.qzv
 
#run PERMANOVA
qiime diversity beta-group-significance \
    --i-distance-matrix ischia_shoot_distance.qza \
    --m-metadata-file shoot_samples_metadata.txt \
    --m-metadata-column Site \
    --p-method permanova \
    --o-visualization Ischia_Site_Significance.qzv
 
#output
#Saved Visualization to: Ischia_Site_Significance.qzv

####converting biom to qza:
(qiime2-2020.2) Winnis-MacBook-Pro:qiime_stuff winni$ qiime tools import --input-path KO_biom.biom --type 'FeatureTable[Frequency]' --output-path KO_feature_table.qza --input-format BIOMV100Format
#instead of --source-format, it’s --input-format
