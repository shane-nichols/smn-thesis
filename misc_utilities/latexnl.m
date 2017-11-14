function stringOut = latexnl(symExprs)
% latex new line (latexnl)

% takes the symmbolic expression symExps and returns the LaTex
% representation. Unix linefeeds '\n' are automatically inserted after
% '\\'.
% \\ is a linefeed in LaTex syntax, but to display them in the Matlab
% Command Window, we must put in \n

lay = latex(symExprs);
lay=strrep(lay,'\','\\');
lay=strrep(lay,'\\\\','\\\\ \n');
stringOut = sprintf(lay);
end