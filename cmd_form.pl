#!/usr/local/bin/perl

use CGI;

print "Content-type: text/html\n\n";

print qq!
<html>
<head>
</head>
<body>
<form action="cmd.pl">
<textarea rows=4 cols=80 name=cmd></textarea>
<br>
<input type=submit>
</form>
</body>
</html>
!;

