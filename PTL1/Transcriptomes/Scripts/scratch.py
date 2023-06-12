import re


og = 'OG6_123456'

ogv = og.split(re.split('OG.{1}_', og)[1])[0][-4:]
print(ogv, re.split('OG.{1}_', og)[1])