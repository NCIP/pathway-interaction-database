#!/usr/local/bin/perl

$cmd = "echo line 1 > carl.out 2>&1";
system($cmd);
$cmd = "echo line 2 >> carl.out 2>&1";
system($cmd);
$cmd = "echo line 3 >> carl.out 2>&1";
system($cmd);
