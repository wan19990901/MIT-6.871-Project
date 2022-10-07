# *Prostate NLP extraction*
# Still Developed, will Add documentation on how to use these code and how the rules are developed

# *Abstract*

## *Objective*

Grade group is considered a major variable of interest when building the a model on prostate cancer. While it's often hard to direcly obtain the grade group or gleason score, the development of natural language proecssing and the abduance of Prostate biopsy pathology report data make it possible for us to extract information without too much manual work. We proposed a program here to to accurately patients' grade group (or gleason score)using information from their biopsy pathology report.

## *Data source*

Research Patient Data Registry (RPDR) Centralized clinical data registry, or data warehouse, that gathers clinical information from various Mass General Brigham hospital systems Queried in September 2020 for patient with diagnosis code or reason for visit of prostate cancer at any time. Query retrieved 88,212 unique patients.

Link : https://www.dropbox.com/s/w5m6qm4wous7ao6/df_pathology_biopsy_final.csv?dl=0

## *Materials and Methods*

A rule based method was developed in the python program in order to extract the grade group from the text. Specifically, the rule is as follow:

### *1: Reports Specification*

We did a filter on report at the very first step. Specifically, we exclude all reports containing the phrase "prostatectomy" or "bone marrow" and then we include the report that contains the both the substring "prost" and one of the substring from "needle","core", or "biops". Also the report should be surgical pyhology report as indicated in the dataset description.

(*Note that this rule is deleted is the newest version as it seems to only appy to the dataset intae was using*)

### *2: Grade Group and Gleason Score*

The definition is that the overall grade group equals to maximum grade group of the individual cores as the Clinical decisions are primarily based on the maximum.

The grade of the core may be described as:
Individual sums of the primary grade (x) and secondary grade (y)
When the primary grade (x) and secondary grade (y) are reported, a grade group is assigned as 1, if Gleason Score 6 (or 3 + 3 = 6).A grade group is assigned as 2, if Gleason Score 7 (or 4 + 3 = 7).A grade group is assigned as 3, if Gleason Score 7 (or 3 + 4 = 7).A grade group is assigned as 4, if Gleason Score is 8 (or 4 + 4 = 8). The other situation where gleason score is 9 or 10 will be considered group 5.

### *3: Extraction Rule:*

#### *Data Preprocesing:*

get all the reports with the following keywords : needle, core, biops, prost,exclude all the reports with the following keywords : prostatectomy, bone marrow, as these two steps would give us all the reports with biopsy.

Find the word gleason; we need to exclude false positive triggered # by M.D. Gleason and found directional word like left or right(Left, RIght, L,R)
\
Once the condition is met: 

get context around word "gleason" (heuristically choose 15). Gleason score we're after should be included in the next 15 words（might need to change this)

remove text in () or []: e.g. Score 3 (60%) + 4 (40%) = 7/10 -> Score 3 + 4 = 7/10,
check if there's only one slash without any + or - signs. 


After after all of the above steps,
we will be able to detect 7 different ways to represent gleason scores. See below :



* 1: when overall grade is represented as primary + secondary

* 2: when overall grade is represented as Gleason 3-4/5

* 3: when overall grade is represented as Gleason 6/10 

* 4: when it is represented as Gleason grade 3/5 and 2/5 (where 3 is the primary and 2 is the secondary grade)

* 5: when it is represented as gleason 3/5 (where primary 3 and secondary 3)confined to..." in the context

* 6: when it is represented as gleason 3 of 5 (where primary 3 and secondary 5) or 6 out of 10 (overall just 6)" in the context

* 7: when it is represented as z (x + y), where z is the total gleason score, x is the primary grade, and y is the secondary grade

Sample examples for the real data:

gleason score x+y  \
gleason score x + y  \
gleason grade x + y  \
gleason x + y  \
gleason 3-4/5 \
gleason grade 3-4/5 \ 
gleason grade 4-5/5 \

gleason score of z \

gleason grade z prostatic adenocarcinoma(not implemented) \
score (x+y)=z (such as score (3+4) = 7)

gleason grades 3/5 and 4/5 \
gleason z/10 \ 
gleason 3/5 \
(gleason grade 3/5) \ 
gleason grade 3/5, \
adenocarcinoma, gleason grade 3/5, present \

Otherwise: 

(Note that the Report_text would be the column for the report)

If Report_text only contains ‘well-differentiated’ (without mention of grade), assign overall grade group = 1 \

If Report_text only contains ‘ moderately differentiated’ (without mention of grade), assign overall grade group = 2
\
If Report_text only contains ‘moderately to poorly differentiated’ (without mention of grade), assign overall grade group = 3
\
If Report_text only contains ‘poorly differentiated’ (without mention of grade), assign overall grade group = 4
\

If text contains ‘gleason score is not considered reliable after hormonal treatment and thus not rendered’, assign overall_grade = 97


## *Instructions*

First Clone the code using git clone command

### *Run directly With Python Scripts(Still in Development)*

### *Run interactively on Python Jupyter Nobteook*

1: Create a notebook on the same directory as the biopsy_grade_extractor.py file
2: Read the report file using pandas, note that the report file should have a column 
named 'Report_Text' as the column for the pathology report
3: import the extractor using from biopsy_grade_extractor import BiopsyGradeExtractor
4: Fit extractor with the data

See more details and examples on the test.ipynb code

## *Experiments and Results*

To be updated by Madhur @ https://github.com/madhurnayan

## *Addtional Materials:*

Main doc : https://docs.google.com/document/d/1yJQBFc7R2qA8OsyvmL7ujQLOSHl3u4C9hyqrJxGTlhU/edit

PSA abstraction from progress notes : https://docs.google.com/document/d/1byko7GQjivx3K-Z-RBkr4frgT8l9N0T0d7SO0wFXKA0/edit?usp=sharing

Biopy Grade Extrators: https://docs.google.com/document/d/1aUFjvz8bumhCnSUDJ8w4CtKoNmGBKYTGeq0cMqoOOl0/edit
