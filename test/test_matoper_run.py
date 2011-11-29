#!/usr/bin/python

import subprocess
import sys
import getopt


def run(bin):
    proc = subprocess.Popen(bin,stdout=subprocess.PIPE)
    s = proc.communicate()[0]
    return float(s)


def compile(compiler,opt,optim,bin):
    args = [compiler,'-o',bin,opt,'-' + optim,'../matoper.F90','test_matoper.F90','-llapack','-lblas']
    args = [a for a in args if len(a) > 0]
    print ' '.join(args)
    subprocess.call(args)



def htmltable(headers,rows,table):
    html = []
    html += ['<table>']

    html += ['<tr>']
    for h in headers:
        html += ['<th>' + h + '</th>']
    html += ['</tr>']


    for cells,row in zip(table,rows):
        html += ['<tr>']
        html += ['<td>' + row + '</td>']

        for cell in cells:
            html += ['<td>' + str(cell) + '</td>']

        html += ['</tr>']

    html += ['</table>']

    return '\n'.join(html)
    


def benchmark(output=None,comp=False):
    table = []

    compilers = ['ifort','pgf90','gfortran']
    optimizations = ['O0','O2','O3']
    options = ["","-DWITH_OPERATORS"]


    headers = ['Compiler']

    for optim in optimizations:
        for opt in options:
            headers.append(optim + ' ' + opt)


    for compiler in compilers:
        row = []

        for optim in optimizations:
            for opt in options:
                bin="./benchmark_matoper_" + compiler + "_" + opt + "_" + optim 
                #print "bin ",bin

                if comp:
                    compile(compiler,opt,optim,bin)
                else:
                    time = [run(bin) for i in range(5)]
                    time = sorted(time)
                    
                    # min max and median
                    print bin, ':' , time[0],time[-1],time[len(time)/2]

                    row.append('%.3f' % time[len(time)/2])


        table.append(row)


    if not comp:
        out = htmltable(headers,compilers,table)

        if output:
            f = open(output,'w')
            f.write(out)
            f.close()
        else:
            print out
    

def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "co:", ["compile", "output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    output = None
    comp = False
    
    for o, a in opts:
        if o in ("-c","--compile"):
            comp = True
        elif o in ("-o", "--output"):
            output = a
        else:
            assert False, "unhandled option"


    benchmark(output,comp)
            
        

if __name__ == '__main__':
    main()
