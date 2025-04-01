# ESMfold_hallucinations

### Project Roadmap

Overview: Identify why protein language models hallucinate, by generating a series of folds that are hallucinations, then evaluating them. Model should be scalable to be able to handle largeamounts of generated folds and initial folds.


	

1. Choose a protein language model an try to make it hallucinate
	
 a. Get Pure language model - ESMfold2 running on computer, HPC cluster etc 
	
 b. Pick three folds (common ones such as proteases, globulins, kinases etc), with maximum of 50% similarity. Pick 10 per category (Use Interpro to export a fasta of proteins, and 	cluster them, choosing far distances) (0-2 resolution filter)
	
 c. Methods to induce changes in existing structures: Maximum likelihood per residue, Temperature, Maximum sequence likelihood (choosing changes that increase likelihood of sequence)

 d. Utilise Markov model or tree to keep track of changes 
	
 e. Generate hallucinations, check against known structures to see if they're hallucinations.

2.	Evaluating hallucinations vs real folds


### Currently incomplete, next steps to be determined based on results of previous step.
