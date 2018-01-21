function [ ] = export_img_latex( figure_handle, filename )

print(figure_handle, filename, '-dpng', '-r0')

end

