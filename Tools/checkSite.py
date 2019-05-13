
import os


with open("fileAODDY2018.txt") as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
content = [x.strip() for x in content] 

print content[0]
b=0
for a in content:
      b=b+1
      command='dasgoclient --query="parent dataset=%s site=T2_US_Purdue" > raw_%s'%(a,b)
      #os.system(command) 
      print command
