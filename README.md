# ESMfold_hallucinations
# This readme is currently incomplete, but it provides an overview and general structure for project, as this is an MSc thesis project in early stages.
### Project Roadmap

Overview: Identify why protein language models hallucinate, by generating a series of folds that are hallucinations, then evaluating them. Model should be scalable to be able to handle largeamounts of generated folds and initial folds.


The general pipeline for this project will be as follows:	

1. Select three protein families, and ten sequences of each of these. These must be as distant as possible.
2. Identification of the maximally distant sequences - this is performed by an embedding and clustering in this pipeline.
3. Generate protein sequences using a tool such as EvoProtGrad with parameters designed to produce hallucinated proteins.
4a. Use sequences from previous pipelines as inputs for structure from sequence models.
4b. ESM3 can also be an option, as it takes prompts from other sources than sequence.
4c. Exploratory data analysis on wild-type/final variant pairs to identify features the model optimises for
5. Use Foldseek to compare generated sequences against known sequences to determine if they are hallucinated.
6. Development of evaluation metrics to identify hallucinations.


