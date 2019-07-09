# Visualize structural variation in a genome

This example code visualizes structural variation in a genome using the BioPython library GenomeDiagram.

We will compare two plasmids (PlasmidA.gbk and PlasmidB.gbk) to find which regions are shared between them.

First, use minimap2 to identify regions that are shared between the two plasmids:

`minimap2 ./PlasmidA.fna ./PlasmidB.fna > PlasmidA_vs_PlasmidB.minimap2.tab`

then filter the output to only include shared regions greater than 10kb for viewing:

`./filter_short.sh`

Next, combine the annotated genbank files for both plasmids:

`cat PlasmidA.gbk PlasmidB.gbk > combined.gbk`

The file "annotationtable.txt" contains annotation that you would like colored on each of the plasmids, and colors are loaded from "colors.txt"


Make the final figure by running the following command:

`python ./run_genomediagram.py ./combined.gbk ./annotationtable.txt --colorfile ./colors.txt --size mbio_half_height --blast PlasmidA_vs_PlasmidB.minimap2.tab.10kb.tab --outname PlasmidA_vs_PlasmidB.pdf --outtype PDF`

The final output in pdf form will show the map:

