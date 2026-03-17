
paper.pdf: paper.md paper.bib
	apptainer run --env JOURNAL=joss docker://openjournals/inara

clean:
	rm -rf paper.pdf jats
