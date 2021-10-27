import React, { useState, useEffect } from 'react';

const BlastDB = (props) => {

    const [data, setData] = useState(null)

    useEffect(()=>{
      console.log('[BlastDB] useEffect');
      setTimeout(()=>{
        setData([{'name': 'nt'}, {'name': 'nr'}])
      }, 6000)
    }, [])
  
    console.log('[BlastDB]', props, data);
    let body = null;
    if (data){
      body = data.map((row)=><div> {row.name}</div>)
    }
    
    if (!body) {
      return <h2>JAWs Blast Database Inventory:</h2>;
    } else {
      return body;
    }
}

// class BlastDB extends Component {
//   state = { data : null }

//   componentDidMount() {
//     console.log('[BlastDB] Ready');
//     setTimeout(()=>{
//       this.setState({data: [{'name': 'nt'}, {'name': 'nr'}]})
//     }, 3000)
//   }
  
//   render() { 
//     console.log('[BlastDB]', this.props, this.state);
//     let body = null;
//     if (this.state.data){
//       body = this.state.data.map((row)=><div> {row.name}</div>)
//     }
    
//     if (!body) {
//       return <h2>Header</h2>;
//     } else {
//       return body;
//     }
//   }
// }
 
export default BlastDB;