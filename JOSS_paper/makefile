
paper.pdf: paper.md
	apptainer run --env JOURNAL=joss docker://openjournals/inara

clean:
	rm -rf paper.pdf jats
