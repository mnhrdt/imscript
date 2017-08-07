cat ../i.html |\
	awk '/<pre class="script">/{b=1;next} /<\/pre>/{b=0} b {print}' |\
	parallel
