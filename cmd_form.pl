#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


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

