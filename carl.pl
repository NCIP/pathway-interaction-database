#!/usr/local/bin/perl


# Copyright SRA International
#
# Distributed under the OSI-approved BSD 3-Clause License.
# See http://ncip.github.com/pathway-interaction-database/LICENSE.txt for details.


$cmd = "echo line 1 > carl.out 2>&1";
system($cmd);
$cmd = "echo line 2 >> carl.out 2>&1";
system($cmd);
$cmd = "echo line 3 >> carl.out 2>&1";
system($cmd);
