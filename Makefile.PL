use inc::Module::Install;

name 'Coursera-NLP-HMM';
all_from 'lib/Coursera/NLP/HMM.pm';

requires 'Data::Dump';
requires 'FindBin';

auto_install;

WriteAll;

