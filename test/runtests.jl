using OncoDiversity
using Test
using DataFrames
import DataFrames: rename!
import OncoDiversity: CDR3SeqAnalysis, countframe, rename!
using CSV
using Query

testcsv(path::AbstractString; kwargs...) = CSV.read(path, DataFrame; kwargs...)


@testset "OncoDiversity.jl" begin
    # Write your tests here.
    @testset "CDRAnalysis" begin
    A = CDR3SeqAnalysis()
    zoodf = testcsv("zoo.csv", header=1)
    @test size(zoodf, 1) == 34
    rename!(zoodf, "Filename", :patient)
    zoocounts = countframe(A, zoodf)
    pd = patientdiversity(zoocounts, 100)
    @test pd.patient == ["zoo1","zoo2"]
    @test pd.diversity == [4.312571345903706, 2.141241195861562]
    @test pd.numclust == [6,5]
    @test pd.checkrange == [1,1]

    pd = patientdiversity(zoocounts, 0.01)
    @test pd.patient == ["zoo1","zoo2"]
    @test pd.diversity == [5.998197562121245, 4.988521617150489]
    @test pd.numclust == [6,5]
    @test pd.checkrange == [1,1]

    #joinkey(i) = (i.patient, i.CDR3)
    # bigframe = @from i in zoocounts begin
    #     @join j in zoodf on i.CDR3 equals j.CDR3
    #     @collect DateFrame
    # end
    # @show bigframe
    @test innerjoin(zoocounts, zoodf, on=[:patient, :CDR3]) |> size == (34,5)
    ds = CDR3SeqData(zoodf)
    @test DataFrame(CDR3SeqData(zoodf)) |> size == (11,5)
    @test summarize(ds) |> size == (2,7)

    qrange =exp10.(range(-2.0, stop=2.0, length=10000)) 
    divs = DiversityScores(zoocounts, qrange)
    @show divs.data
    @test 7 < divs.data.IPq[1] < 8
    @test 1.5 < divs.data.IPq[2] < 1.9
    @test 0.5 < divs.data.IPslope[1] < 0.6
    @test 1.0 < divs.data.IPslope[2] < 1.05
    @test 1.7 < divs.data.Shan[1] < 1.8
    @test 1.3 < divs.data.Shan[2] < 1.4
    @test 0.1 < divs.data.Simp[1] < 0.2
    @test 0.25 < divs.data.Simp[2] < 0.35
    
    ## adding in a "skewed zoo"
    zoodf2 = testcsv("zoo2.csv", header=1)
    @test size(zoodf2, 1) == 17*3
    rename!(zoodf2, "Filename", :patient)
    zoocounts2 = countframe(A, zoodf2)
    pd2 = patientdiversity(zoocounts2, 100)
    pd2 = patientdiversity(zoocounts2, 0.01)

    @test innerjoin(zoocounts2, zoodf2, on=[:patient, :CDR3]) |> size == (17*3,5)
    ds2 = CDR3SeqData(zoodf2)
    @test summarize(ds2) |> size == (3,7)

    divs2 = DiversityScores(zoocounts2, qrange)
    @show divs2.data

    ## mimicking low recovery dataset
    zoodf3 = testcsv("zoo3.csv", header=1)
    zoocounts3 = countframe(A, zoodf3)
    pd3 = patientdiversity(zoocounts3, 100)
    pd3 = patientdiversity(zoocounts3, 0.01)
    ds3 = CDR3SeqData(zoodf3)
    divs3 = DiversityScores(zoocounts3, qrange)
    @show divs3.data

end
end
