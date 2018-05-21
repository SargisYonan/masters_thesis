function [ ] = export_img_latex( figure_handle, filename )

type = 'png';
dpi = '1000';

filename = strcat([filename, '.', type]);
fileext = strcat(['-d', type]);

res = strcat(['-r', dpi]);

print(figure_handle, filename, fileext, res);

end

