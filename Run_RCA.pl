
#!/usr/bin/perl 

# my $dir= "C:\Documents and Settings\yanbi.NIH\Desktop\matlab\RCA";
# chdir($dir);
# $array = "'input_array.txt'";
# $TF = "'5TF_lncR_input_binding.txt'";

$array = "'mRNA13_up.txt'";
$a = "'median'";

$TF = "'5mir_input_binding.txt'";
$cmd = 'matlab -nodisplay -wait -nosplash -nodesktop -r "RCA('.$array .', '.$TF.', '.$a.', 0, 1);quit;"';
system($cmd);