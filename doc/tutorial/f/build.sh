cat ../i.html | awk '/<\/pre>/{x=0}x;/<pre class="script">/{x=1}' | parallel -j 5
