#reset
#set term pdf
#set out 'fig4.pdf'

set grid

plot 'fort.10' w l t'Interp.', cos(x) t'cos(x)', 'fort.11' w p t'Points'
