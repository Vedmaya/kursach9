function data = readFile( filename, x, y ) %x-strings, y-columns
    file=fopen(['backup/' filename],'r');
    data=fread(file,[x,y],'float64');
    fclose(file);
end

