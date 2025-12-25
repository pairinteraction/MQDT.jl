using MQDT
using PythonCall

function test_model_struct(T::fModel)
    if (T.size != length(T.terms)) || (T.size != length(T.core)) || (T.size != size(T.defects, 1))
        println("Model size does not correspond to provided parameters.")
        @test false
    elseif size(T.unitary) != (T.size, T.size)
        println("Frame transformation matrix does not have the correct dimensions.")
        @test false
    elseif length(findall(T.core)) != size(T.outer_channels) != size(T.inner_channels)
        println("Wrong number of channels.")
        @test false
    elseif length(T.mixing) != size(T.angles, 1)
        println("Wrong number of close-coupling rotations.")
        @test false
    else
        println("Model $(T.name) passed.")
        @test true
    end
end

function test_model_struct(T::kModel)
    if (T.size != length(T.terms)) ||
       (T.size != length(T.jjscheme)) ||
       (T.size != length(T.lschannels.i)) ||
       (T.size != length(T.jjchannels.i)) ||
       (T.size != length(T.K1))
        println("Model size does not correspond to provided parameters.")
        @test false
    elseif size(T.K0) != (T.size, T.size)
        println("Frame transformation matrix does not have the correct dimensions.")
        @test false
    else
        println("Model $(T.name) passed.")
        @test true
    end
end

function test_model_struct(T::Union{Vector{fModel},Vector{kModel}})
    for i in eachindex(T)
        test_model_struct(T[i])
    end
end

function test_model_unitary(T::fModel, P::Parameters)
    c = findall(T.core)
    t0 = T.unitary[c, c]

    if any(qn -> isnan(qn.Jc), T.outer_channels.i) ||
       (typeof(T.outer_channels) == fjChannels && any(qn -> isnan(qn.Fc), T.outer_channels.i))
        println("Skip test_model_unitary for model $(T.name) due to NaN in quantum numbers.")
        return
    end

    if typeof(T.outer_channels) == jjChannels
        if typeof(T.inner_channels) == lsChannels
            t1 = MQDT.matrix_ls_to_jj(T.inner_channels, T.outer_channels)'
        else
            s = size(T.inner_channels)
            t1 = diagm(repeat([1.0], s))
        end
    elseif typeof(T.outer_channels) == fjChannels
        if typeof(T.inner_channels) == jjChannels
            t1 = MQDT.matrix_jj_to_fj(T.inner_channels, T.outer_channels, P.spin)'
        else
            lc = unique(MQDT.get_lc(T.inner_channels))
            lr = unique(MQDT.get_lr(T.inner_channels))
            ft = unique(MQDT.get_F(T.outer_channels))[1]
            qn = MQDT.AngularMomenta(lc, lr, P.spin)
            jj = MQDT.jj_channels(qn)
            fj = [unique.(MQDT.get_F.(jj))[i][1] for i in 1:length(jj)]
            i = findfirst(isequal(ft), fj)
            tjj = MQDT.matrix_ls_to_jj(T.inner_channels, jj[i])
            tfj = MQDT.matrix_jj_to_fj(jj[i], T.outer_channels, P.spin)
            t1 = (tjj * tfj)'
        end
    end

    @test isapprox(t0, t1; rtol=0.001)
end

function test_model_unitary(T::kModel)
    # TODO not working yet due to NaNs in quantum numbers
    println("Skip test_model_unitary for kModel $(T.name).")
    return
    t1 = eigen(T.K0).vectors
    t2 = MQDT.matrix_ls_to_jj(T.lschannels, T.jjchannels)'
    @test isapprox(t1, t2; rtol=0.001)
end

function test_model_states_and_matrix_elements(model::fModel, param::Parameters)
    nu_min, nu_max = get_nu_limits_from_model(model)
    nu_max = isnan(nu_max) ? nu_min + 5 : min(nu_max, nu_min + 5)
    states = eigenstates(nu_min, nu_max, model, param)
    basis = basisarray(states, model)
    @test size(basis) == length(states)

    logging = pyimport("logging")
    logging.getLogger("").setLevel(logging.ERROR)
    me = matrix_elements(basis, param)
    @test length(me["dipole"]) == 0
    @test length(me["quadrupole"]) > 0 || length(me["diamagnetic"]) > 0

    qns0 = model.outer_channels.i[1]
    if qns0.F != 0
        @test length(me["paramagnetic"]) >= size(basis)
    else
        @test length(me["paramagnetic"]) == 0
    end
end

@testset "Models" begin
    for species in (MQDT.Sr87, MQDT.Sr88, MQDT.Yb171, MQDT.Yb173, MQDT.Yb174)
        for nm in names(species, all=true)
            if !isdefined(species, nm)
                continue
            end
            obj = getfield(species, nm)
            if !isa(obj, MQDT.Model)
                continue
            end

            println("Testing $(species) $(nm)")
            @testset "$(species) $(nm)" begin
                if isa(obj, fModel)
                    test_model_struct(obj)
                    test_model_unitary(obj, species.PARA)
                    test_model_states_and_matrix_elements(obj, species.PARA)
                elseif isa(obj, kModel)
                    test_model_struct(obj)
                    test_model_unitary(obj)
                end
            end
        end
    end
end
