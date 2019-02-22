
#!/usr/bin/perl 

# my $dir= "C:\Documents and Settings\yanbi.NIH\Desktop\matlab\RCA";
# chdir($dir);

# $array = "'input_array.txt'";
$array = "'down.txt'";
$a = "'median'";

# $TF = "'input_TF_binding.txt'";
$TF = "'new_input_TF_binding.txt'";
$cmd = 'matlab -nodisplay -wait -nosplash -nodesktop -r "RCA('.$array .', '.$TF.', '.$a.', 0, 1);quit;"';
system($cmd);