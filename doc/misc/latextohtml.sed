s_<h1>\(.*\)</h1>_{\\LARGE \1}_g
s_<h2>\(.*\)</h2>_\\section{\1}_g
s_<h3>\(.*\)</h3>_\\subsection{\1}_g
s_<h4>\(.*\)</h4>_\\subsubsection{\1}_g
s_<p>__g
s_</p>__g
s_<i>\([^<]*\)</i>_\\emph{\1}_g
s_<b>\([^<]*\)</b>_{\\bf \1}_g
s_<i>_\\emph{_
s_</i>_}_
s_<b>_{\\bf _
s_</b>_}_
s_<table.*>_\\begin{tabular}{l|l}_
s_</table>_\\end{tabular}_
s_</td>.*</tr>_\\\\_g
s_</th>.*</tr>_\\\\\\hline_g
s_</td>_\&_g
s_</th>_\&_g
s_<tr>__g
s_<td>__g
s_<th>__g
s_<br.*>__g
s_<blockquote>_\\begin{quotation}_
s_</blockquote>_\\end{quotation}_
s_{% highlight \(.*\) %}_\\begin{minted}{\1}_
s_{% endhighlight %}_\\end{minted}_
s_section{[0-9.]*\s_section{_
