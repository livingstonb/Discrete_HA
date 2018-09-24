function variance = variance_discrete(probabilities,values)

Evalues = probabilities' * values;
variance = probabilities' * (values - Evalues).^2;

end