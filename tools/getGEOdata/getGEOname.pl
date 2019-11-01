#perl parser to get title and file size from GEO accession number

$geo = $ARGV[0];

my $url = `curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$geo"`;

if( $url =~ /<td nowrap>Title<\/td>\s*<td\sstyle="text-align:\sjustify">(.+)<\/td>/) {
	print STDOUT $1;

	# check for file size
	if( $url =~ /<td><a\shref="ftp(.+)"\starget="_blank">SOFT\sformatted\sfamily\sfile/) {

		$dataURL = `curl -s https$1`;
		if( $dataURL =~ /\s+(\d+\.?\d*[M,K,G])\s+<hr>/ ) {
			print STDOUT " [data size: $1]";
		}
	}

} else {
	print STDOUT "no match";
}





