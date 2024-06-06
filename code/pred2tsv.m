% "convert prediction in the MATLAB data structure pred to a tsv"
% Usage: pred2tsv(pred_path, outputFileName)

function pred2tsv(pred_path, outputFileName)
    load(pred_path)
    % Open the file for writing
    fileID = fopen(outputFileName, 'w');

    for p = 1:size(pred.object, 1)
        % Access each protein
        current_protein = pred.object{p};

        % Find non-zero scores for the current object
        non_zero_indices = full(pred.score(p, :)) ~= 0;

        % Preallocate a structure array for predicted terms
        predicted_terms = struct('term', {}, 'score', {});

        % Extract terms and scores in a structure array
        terms = pred.ontology.term(non_zero_indices);
        scores = full(pred.score(p, non_zero_indices));

        for idx = 1:length(terms)
            predicted_terms(idx).term = terms(idx);
            predicted_terms(idx).score = scores(idx);
        end

        % Write the predicted terms to the file
        for idx = 1:length(predicted_terms)
            term_value = predicted_terms(idx).term.id; % Extract the GO term
            score_value = predicted_terms(idx).score; % Extract the score

            % Print the current object, term, and score to the file
            fprintf(fileID, '%s\t%s\t%.10f\n', current_protein, term_value, score_value);
        end
    end

    % Close the file
    fclose(fileID);

    disp(['Predicted terms written to ', outputFileName]);
end

