import re, sys


def readVowelsOnlyDict():
    vowfile = open('phondicts/vowels_only.dat', 'r')
    vowfile.readline()
    vdict = {}
    for line in vowfile:
        re.split(',', line)
    
    
