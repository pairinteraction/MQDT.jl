using MQDT


function test_model(T::fModel)
    if T.size != length(T.terms) != length(T.core) != length(T.thresholds) != size(T.defects, 1)
        return println("Model size does not correspond to provided parameters.")
    elseif size(T.unitary) != (T.size, T.size)
        return println("Frame transformation matrix does not have the correct dimensions.")
    elseif length(findall(T.core)) != size(T.outer_channels) != size(T.inner_channels)
        return println("Wrong number of channels.")
    elseif length(T.mixing) != size(T.angles, 1)
        return println("Wrong number of close-coupling rotations.")
    else
        return println("Model $(T.name) passed.")
    end
end

function test_model(T::kModel)
    if T.size !=
       length(T.terms) !=
       length(jjscheme) !=
       length(T.lschannels) !=
       length(T.jjchannels) !=
       length(T.thresholds) !=
       length(T.K1)
        return println("Model size does not correspond to provided parameters.")
    elseif size(T.K0) != (T.size, T.size)
        return println("Frame transformation matrix does not have the correct dimensions.")
    else
        return println("Model $(T.name) passed.")
    end
end

function test_model(T::Union{Vector{fModel},Vector{kModel}})
    for i in eachindex(T)
        test_model(T[i])
    end
end

function test_unitary(T::fModel, P::Parameters)
    c = findall(T.core)
    t0 = T.unitary[c, c]
    if typeof(T.outer_channels) == jjChannels
        if typeof(T.inner_channels) == lsChannels
            t1 = matrix_ls_to_jj(T.inner_channels, T.outer_channels)'
        else
            s = size(T.inner_channels)
            t1 = diagm(repeat([1.0], s))
        end
    elseif typeof(T.outer_channels) == fjChannels
        if typeof(T.inner_channels) == jjChannels
            t1 = matrix_jj_to_fj(T.inner_channels, T.outer_channels, P.spin)'
        else
            lc = unique(get_lc(T.inner_channels))
            lr = unique(get_lr(T.inner_channels))
            ft = unique(get_F(T.outer_channels))[1]
            qn = mqdt.AngularMomenta(lc, lr, P.spin)
            jj = jj_channels(qn)
            fj = [unique.(get_F.(jj))[i][1] for i in 1:length(jj)]
            i = findfirst(isequal(ft), fj)
            tjj = matrix_ls_to_jj(T.inner_channels, jj[i])
            tfj = matrix_jj_to_fj(jj[i], T.outer_channels, P.spin)
            t1 = (tjj * tfj)'
        end
    end
    return isapprox(t0, t1), t0, t1
end

function test_unitary(T::kModel)
    t1 = eigen(T.K0).vectors
    t2 = matrix_ls_to_jj(T.lschannels, T.jjchannels)'
    return t1, t2
end