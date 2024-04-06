window.BENCHMARK_DATA = {
  "lastUpdate": 1712415758809,
  "repoUrl": "https://github.com/SciFracX/FractionalCalculus.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "committer": {
            "email": "2283984853@qq.com",
            "name": "ErikQQY",
            "username": "ErikQQY"
          },
          "distinct": true,
          "id": "82a0b46ddcf5724f55d00dc79b77f2be5998371c",
          "message": "Fix benchmarks\n\nSigned-off-by: ErikQQY <2283984853@qq.com>",
          "timestamp": "2024-04-06T23:00:01+08:00",
          "tree_id": "af9a292da9b481d3d361e76139064e93296c021b",
          "url": "https://github.com/SciFracX/FractionalCalculus.jl/commit/82a0b46ddcf5724f55d00dc79b77f2be5998371c"
        },
        "date": 1712415758154,
        "tool": "julia",
        "benches": [
          {
            "name": "Caputo/CaputoL1",
            "value": 4291,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":7,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Caputo/CaputoTrap",
            "value": 7479,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":4,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Caputo/CaputoDiethelm",
            "value": 9178,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}