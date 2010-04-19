import re, sys

vowels = ['AO', 'AA', 'IY', 'UW', 'EH', 'IH', 'UH', 'AH', 'AE', 'EY', 'AY', 'OW', 'AW', 'OY']

def destress(p):
    return re.findall(r'[^012]+', p)[0]

lexdict = {}

infile = open('IPhOD2_Words.txt', 'r')
wordsOutfile = sys.stdout #open('vowels_only.dat', 'w')

lines = re.split(r'\n', infile.read())

outfile.write(','.join(['word', 'phones', 'vowels', 'freq']) + '\n')

for line in lines[1:-1]:
    # split the line into fields
    lsp = re.split(r'\s+', line)
    # split the phonetic spelling into phones
    phons = re.split(r'\.', lsp[3])
    # get the vowels
    vs = [p for p in phons if destress(p) in vowels]
    key = tuple([destress(p).lower() for p in vs])
    try: lexdict[key] += lsp[-3]
    except KeyError: lexdict[key] = lsp[-3]
    wordsOutfile.write(','.join((lsp[1], lsp[3], '.'.join(vs), lsp[-3])) + '\n')

