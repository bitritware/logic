"""Builds the unary table (LaTeX) and the CSV from unary_catalogue.json (written by chai_checker.py).
Usage: python3 make_tables.py [catalogue.json] [out_dir]"""
import json, csv, re, sys, os
cat_path=sys.argv[1] if len(sys.argv)>1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),'unary_catalogue.json')
out=sys.argv[2] if len(sys.argv)>2 else os.path.dirname(os.path.abspath(__file__))
cat=json.load(open(cat_path))
dig={0:'0',1:'1',-1:'2'}; sym={-1:'-',0:'0',1:'+'}
rows=[]
for v,(s7,t7),(s2,t2) in cat:
    rows.append((s7,''.join(dig[x] for x in v),''.join(sym[x] for x in v),t7,s2,t2))
rows.sort(key=lambda r:(r[0],r[1]))
with open(os.path.join(out,'unary_catalogue.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['code_0t','vector','size_seven_ops','term_seven_ops','size_DIF_NOT','term_DIF_NOT'])
    for s7,c,vec,t7,s2,t2 in rows: w.writerow(['0t'+c,vec,s7,t7,s2,t2])
def tex(t):
    t=re.sub(r'\b(SYN|DIF|OR|XOR|AND|NOT|ANT)\b',lambda m:'\\'+m.group(1),t)
    return '$'+t+'$'
def vec(v): return '$'+v.replace('-','-').replace('+','+')+'$'
n=len(rows); half=(n+1)//2
L=rows[:half]; R=rows[half:]
lines=[]
for i in range(half):
    a=L[i]; cells=[r'\texttt{0t%s}'%a[1],vec(a[2]),str(a[0]),tex(a[3])]
    if i<len(R):
        b=R[i]; cells+=[r'\texttt{0t%s}'%b[1],vec(b[2]),str(b[0]),tex(b[3])]
    else: cells+=['','','','']
    lines.append(' & '.join(cells)+r'\\')
body='\n'.join(lines)
open(os.path.join(out,'S5_unary_table.tex'),'w').write(r'''\begin{table*}[t]
\centering\scriptsize
\setlength{\tabcolsep}{3pt}
\caption{The 27 unary functions $\Fone$ (the G3 branches) as \bidra{} terms. Code: digits $0,1,2$ for $0,+1,-1$, listing $g(-1),g(0),g(+1)$, with the prefix \texttt{0t} marking a ternary literal as \texttt{0b} and \texttt{0x} do for binary and hexadecimal. Size: minimum number of operator nodes over the seven operators (tree size, repeated subterms counted each time; exhaustive search). Term: the lexicographically smallest minimal term.}
\label{tab:unary}
\begin{tabular}{llcl@{\hspace{1.2em}}llcl}
\toprule
code & $g$ & size & term & code & $g$ & size & term\\\midrule
'''+body+r'''
\bottomrule
\end{tabular}
\end{table*}
''')
print(n,'rows',half,'left')
