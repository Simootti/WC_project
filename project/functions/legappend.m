function legappend(newStrings)

    % Quick & dirty fork of legappend specific to MATLAB versions >= R2014b
    % Only supports appending strings to the existing legend handle
    % Assumes only one legend is in the current figure
    % Add multiple strings by passing it a 1D cell array of strings

    % Find our legend object
    h = findobj(gcf, 'Type', 'legend');

    if ~isempty(h)
        % Get existing array of legend strings and append our new strings to it
        oldstr = h.String;
        if ischar(newStrings)
            % Input string is a character array, assume it's a single string and
            % dump into a cell
            newStrings = {newStrings};
        end

        newstr = [oldstr newStrings];

        % Get line object handles
        ploth = flipud(get(gca, 'Children'));

        % Update legend with line object handles & new string array
        h.PlotChildren = ploth;
        h.String = newstr;
    end
end