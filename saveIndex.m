function saveIndex( i )
    file=fopen('backup/I','w+');
    fwrite(file,i,'float64');
    fclose(file);
end

