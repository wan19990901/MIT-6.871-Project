# Prostate NLP extraction 

# TODO: Add documentation on how to use these code and how the rules are developed.

# Abstract

## Objective

Grade group is considered a major variable of interest when building the a model on prostate cancer. While it's often hard to direcly obtain the grade group or gleason score, the development of natural language proecssing and the abduance of Prostate biopsy pathology report data make it possible for us to extract information without too much manual work. We proposed a program here to to accurately patients' grade group (or gleason score)using information from their biopsy pathology report.

## Data source
Research Patient Data Registry (RPDR) Centralized clinical data registry, or data warehouse, that gathers clinical information from various Mass General Brigham hospital systems Queried in September 2020 for patient with diagnosis code or reason for visit of prostate cancer at any time. Query retrieved 88,212 unique patients.

## Materials and Methods

A rule based method was developed in the python program in order to extract the grade group from the text. Specifically, the rule is as follow:

### 1: Reports Specification

We did a filter on report at the very first step. Specifically, we exclude all reports containing the phrase "prostatectomy" or "bone marrow" and then we include the report that contains the both the substring "prost" and one of the substring from "needle","core", or "biops". Also the report should be surgical pyhology report as indicated in the dataset description.

### 2: Grade Group and Gleason Score

The definition is that the overall grade group equals to maximum grade group of the individual cores as the Clinical decisions are primarily based on the maximum.

The grade of the core may be described as:
Individual sums of the primary grade (x) and secondary grade (y)
When the primary grade (x) and secondary grade (y) are reported, a grade group is assigned as 1, if Gleason Score 6 (or 3 + 3 = 6).A grade group is assigned as 2, if Gleason Score 7 (or 4 + 3 = 7).A grade group is assigned as 3, if Gleason Score 7 (or 3 + 4 = 7).A grade group is assigned as 4, if Gleason Score is 8 (or 4 + 4 = 8). The other situation where gleason score is 9 or 10 will be considered group 5.

### 3: Extraction Rule:

1: Find patterns of Individual sums of the primary grade (x) and secondary grade (y)

gleason score x+y or
gleason score x + y or 
gleason grade x + y or
gleason x + y or 
gleason grade of x and y(not yet)

2: Total sum (z) of the primary (x) and secondary (y) grades

note: If only total score 7 and no primary or secondary grades available i.e. cannot distinguish grade group 2 or 3, then assign grade group 2.5


gleason score of z
gleason z/10
gleason grade z prostatic adenocarcinoma
score (x+y)=z (such as score (3+4) = 7)

3: Grade group

Grade group 1: Find patterns on :

Gleasong grade (3/5) (Typo is deliberate)
gleason 3/5
gleason score 3 (of 5)
(gleason grade 3/5)
, gleason grade 3/5,
adenocarcinoma, gleason grade 3 out of 5, involving
adenocarcinoma, gleason grade 3/5, present
gleason grade 3 (score 6 of 10)

Grade group 3: Find patterns on :

gleason 3-4/5
gleason grade 3-4/5
gleason grades 3/5 and 4/5

Grade group 4, Find patterns on :

Gleason grade 4/5

Grade group 5, Find patterns on :

gleason grade 4-5/5


Otherwise:

(Note that the Report_text would be the column for the report)

If Report_text only contains ‘adenocarcinoma, well-differentiated’ (without mention of grade), assign overall grade group = 1
If Report_text only contains ‘adenocarcinoma, well to moderately differentiated’ or ‘ moderately differentiated’ (without mention of grade), assign overall grade group = 2
If Report_text only contains ‘adenocarcinoma, moderately to poorly differentiated’ (without mention of grade), assign overall grade group = 3
If Report_text only contains ‘adenocarcinoma, poorly differentiated’ (without mention of grade), assign overall grade group = 4
If text contains ‘gleason score is not considered reliable after hormonal treatment and thus not rendered’, assign overall_grade = 97


## Instructions

First Clone the code using git clone command

### Run directly With Python Scripts(Still in Development)

### Run interactively on Python Jupyter Nobteook

1: Create a notebook on the same directory as the biopsy_grade_extractor.py file
2: Read the report file using pandas, note that the report file should have a column 
named 'Report_Text' as the column for the pathology report
3: import the extractor using from biopsy_grade_extractor import BiopsyGradeExtractor
4: Fit extractor with the data

See more details and examples on the test.ipynb code

## Experiments and Results


## Addtional Materials:

Main doc : https://docs.google.com/document/d/1yJQBFc7R2qA8OsyvmL7ujQLOSHl3u4C9hyqrJxGTlhU/edit

PSA abstraction from progress notes : https://docs.google.com/document/d/1byko7GQjivx3K-Z-RBkr4frgT8l9N0T0d7SO0wFXKA0/edit?usp=sharing
