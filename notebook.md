Proposal:
	The goal for this experiment was to change the workflow of and finish the single and bulk cell simulator for the scsim paper.

Result:
	The experiment (more of a project) was successful and the simulator is simulating bulk samples and single cell samples as well as showing proof of usability by running variant callers such as monovar and bcftools on the simulated samples. The only real trouble I had with the project was snakemakes constraints on wildcards in rules. After writing all of the rules I needed to follow the new workflow, I leared more about what I can do with wildcards and what I cannot, and much of the time I spent on the project was fixing snakemake rules. Also, I ended up using an edited version of Patricks mix-fastq program to make the bulk sample bams, but I wrote my own json generator in the snakemake rule itself, as I didnt think Patricks would be easily changed to do what I needed it to.

Conlclusion:
	The similator now needs to have the hiearchical structure implemented, but that should be an easy last step. Then it will be ready for the scsim paper and the scseq paper next semester.
