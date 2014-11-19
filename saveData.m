function saveData( filename, data )
    file=fopen(['backup/' filename], 'a+');
    fwrite(file,data,'float64');
    fclose(file);
end

