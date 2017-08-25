# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 14:59:11 2016

@author: k
"""

import unicodecsv
folder = "E:\\Study\\coursera\\Udacity\\dataFilesUsedInClass\\"

enrollments_filename = folder+ 'enrollments.csv'

## Longer version of code (replaced with shorter, equivalent version below)

# enrollments = []
# f = open(enrollments_filename, 'rb')
# reader = unicodecsv.DictReader(f)
# for row in reader:
#     enrollments.append(row)
# f.close()

with open(enrollments_filename, 'rb') as f:
    reader = unicodecsv.DictReader(f)
    enrollments = list(reader)
    
### Write code similar to the above to load the engagement
### and submission data. The data is stored in files with
### the given filenames. Then print the first row of each
### table to make sure that your code works. You can use the
### "Test Run" button to see the output of your code.

engagement_filename = folder+ 'daily_engagement.csv'
submissions_filename = folder+ 'project_submissions.csv'

with open(engagement_filename,'rb') as f:
    reader = unicodecsv.DictReader(f)
    daily_engagement = list(reader)
with open(submissions_filename,'rb') as f:
    reader = unicodecsv.DictReader(f)
    project_submissions = list(reader)



len(enrollments)

unique_enrolled_students = set()
for enrollment in enrollments:
    unique_enrolled_students.add(enrollment['account_key'])
len(unique_enrolled_students)

len(daily_engagement)

unique_engagement_students = set()
for engagement_record in daily_engagement:
    unique_engagement_students.add(engagement_record['acct'])
len(unique_engagement_students)

len(project_submissions)

unique_project_submitters = set()
for submission in project_submissions:
    unique_project_submitters.add(submission['account_key'])
len(unique_project_submitters)

count = 0
strange =[]
for enrollment in enrollments:
    student = enrollment['account_key']
    if student not in unique_engagement_students:
        if enrollment['join_date'] != enrollment['cancel_date']:
            count +=1
            strange.append(enrollment)

paid_students = {}
from datetime import ddatetime as dt
for enrollment in enrollments:
    

