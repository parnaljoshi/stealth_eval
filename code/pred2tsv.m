function pred2tsv(pred, outputFileName)
    % Open the file for writing
    fileID = fopen(outputFileName, 'w');

    for p = 1:size(pred.object, 1)
        % Access the current object
        current_object = pred.object{p};

        % Find non-zero scores for the current object
        non_zero_indices = full(pred.score(p, :)) ~= 0;

        % Preallocate a structure array for predicted terms
        predicted_terms = struct('term', {}, 'score', {});

        % Assign terms and scores to the structure array
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
            fprintf(fileID, '%s\t%s\t%.10f\n', current_object, term_value, score_value);
        end
    end

    % Close the file
    fclose(fileID);

    disp(['Predicted terms written to ', outputFileName]);
end

