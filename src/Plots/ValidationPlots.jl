"""
    pp_prob_plot(sepp::SEPP, pp::PP)
Plot the probability plot for the point process.

See 4.1 in Li2020. 
""" 
function pp_prob_plot(sepp::SEPP, pp::PP)
    s, p = pp_analysis(sepp, pp)

    id = layer(x = p, y = p, color = [color("red")], Geom.line, order = 2)

    emp = layer(x = p, y = s, color = [color("black")], Geom.line, order = 1)

    plt = plot(id, emp, Guide.title = "Point process probability plot", Guide.xlabel = "Empirical", Guide.ylabel = "Model")

    return plt
end

"""
    marks_prob_plot(sempp::SEMPPExpKern, mpp::MarkedPointProcess)
Plot the probability plot for the marks.

See 4.2 in Li2020. 
""" 
function marks_prob_plot(sempp::SEMPPExpKern, mpp::MarkedPointProcess) 
    s, p = transformed_marks_ecdf(sempp, mpp)

    id = layer(x = p, y = p, color = [color("red")], Geom.line, order = 2)

    emp = layer(x = p, y = s, color = [color("black")], Geom.point, order = 1)

    plt = plot(id, emp, Guide.title = "Marks probability plot", Guide.xlabel = "Empirical", Guide.ylabel = "Model")

    return plt
end

"""
    marks_qq_plot(sempp::SEMPPExpKern, mpp::MarkedPointProcess)
Plot the quantile plot based on the unit exponential distribution.

See 4.2 in Li2020.
"""
function marks_qq_plot(sempp::SEMPPExpKern, mpp::MarkedPointProcess)
    emp_q, mod_q = marks_unit_exponential_qq(sempp, mpp)

    id = layer(x = mod_q, y = mod_q, color = [color("red")], Geom.line, order = 2)
    
    emp = layer(x = mod_q, y = emp_q, color = [color("black")], Geom.point, order = 1)
    
    plt = plot(id, emp, Guide.title = "Marks quantile plot", Guide.xlabel = "Empirical", Guide.ylabel = "Model")

    return plt
end