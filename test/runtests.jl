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
end
end
