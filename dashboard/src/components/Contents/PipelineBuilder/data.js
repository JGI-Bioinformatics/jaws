const data = {
  "steps": [ 
    { "name": "Start", 
      "level": 0,
      "data": {"status": "completed"}
    },
    { "name": "Filter",
      "level": 1,
      "data": {"status": "completed"}
    },
  //   { "name": "Read QC",
  //     "level": 2,
  //     "branch": { "size": 3, "position": -1},
  //     "data": {"status": "completed"}
  //   },
  //   { "name": "Assembly",
  //     "level": 2,
  //     "branch": { "size": 3, "position": 0},
  //     "data": {
  //       "status": "completed",
  //     }
  //   },
  //   { "name": "Kmer Analysis",
  //     "level": 2,
  //     "branch": { "size": 3, "position": 1},
      
  //     "data": {
  //     "status": "completed",
  //    }
  //   },
  //   {
  //     "name": "Post Analysis",
  //     "level": 3,
  //     "data": {"status": "failed"}
  //   },
  //   {
  //     "name": "Jamo Prep",
  //     "level": 4,
  //     "data": {
  //       "status": ""
  //     }
  //   }
  ],
  "workflow": [
      {"source": "Start", "target": "Filter"},
      // {"source": "Filter", "target": "Assembly"},
      // {"source": "Filter", "target": "Read QC"},
      // {"source": "Filter", "target": "Kmer Analysis"},
      // {"source": "Assembly", "target": "Post Analysis"},
      // {"source": "Read QC", "target": "Post Analysis"},
      // {"source": "Kmer Analysis", "target": "Post Analysis"},
      // {"source": "Post Analysis", "target": "Jamo Prep"}
    ]
}

export default data;
