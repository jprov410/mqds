from keyword_library import keywords,isfloat,default
infile='run.in'
outfile='processed_run.in'

"""
Find the keywords defined in the input file
"""

located={}
with open(infile, 'r') as input:
    for line in input:
        for word in keywords:            
            if word == line.split()[0]:
                if ( line.split()[1].isdigit() ) == True:
                    located.update( { line.split()[0] : int( line.split()[1] ) } )
                elif ( isfloat(line.split()[1]) ) == True:
                    located.update( { line.split()[0] : float( line.split()[1] ) } )
                else:
                    located.update( { line.split()[0] : line.split()[1] } )

"""
If not found in input file, add them with their default valued
that are defined in keyword_library
"""

for word in keywords:
    if word not in located:
        located.update( {word : default(word) })

"""
write the new input in standardized format
"""

with open(outfile,'w') as output:
    for word in sorted(located):
        output.write( word + '\t' + str( located[word]) + '\n' )
