/*
Ref: https://cloudblogs.microsoft.com/opensource/2019/02/28/tutorial-adding-distributed-tracing-instrumentation-to-browser-javascript-app/
span ID: an arbitrary unique ID used to identify this particular unit of work
trace ID: an arbitrary unique ID describing the trace this span is a part of. In a typical backend application, a trace might represent a single user HTTP request to a web application. In the browser, it’s more likely to represent a larger unit of work, like a page view.
parent span ID: this field is used to draw a causal relationship between two spans. For example, if your web application makes an http request to a backend service that then triggers a database query, the span for the database query would use the span ID of the backend service span as its parent span ID.
name: the particular type of work this span represents, e.g. “render json,” “select,” “save to S3,” etc.
service name: the area of the codebase this span came from. In the browser, this might help you identify a particular package or area of your code (“router,” “redux,” “async renderer”).
duration: the amount of time this particular unit of work took, in milliseconds.
*/

import Libhoney from "libhoney";

const honeycomb = new Libhoney({
                      writeKey: process.env.REACT_APP_HONEYCOMB_KEY, // Key from env variable
                      dataset:
                      process.env.NODE_ENV === "production"
                      ? "jaws-dashboard-prod"
                      : "jaws-dashboard-dev"
                    });
 
/** Create a new string ID that is a valid Honeycomb trace/span ID. Currently no gaurentee of uniqueness. */
export const hcCreateId = () => Math.floor(Math.random() * 2 ** 31).toString();

 /**
 * Build up a trace object for honeycomb containing the 3 tracing ID's: trace, span, and parent
 * @param  {string} traceID The ID of the trace that the new span belongs to
 * @param  {string} spanID The ID of the new span, if not provided a new ID is generated
 * @param  {string} parentID The ID of the spans direct parent, can be null for root level spans
 * @returns An object containing all trace ids with honeycomb consistent id casing
 */
export const hcBuildTrace = (traceID, spanID = this.spanIdGen(), parentID) => ({
  trace_id: traceID,
  span_id: spanID,
  parent_id: parentID
});

export const hcRootTraceId = hcCreateId();

/** Send a span with the provided payload and tracing IDs object to the JGI Honeycomb */
export const hcSendSpan = (payload, trace) => {
  const honeycombPayload = {
    // Span Tracing IDs
    // trace: { ...trace },
    "trace.trace_id": trace.trace_id,
    "trace.span_id": trace.span_id,
    "trace.parent_id": trace.parent_id,

    // Send the rest of our fields
    ...payload
  };

  const event = honeycomb.newEvent();
  // console.log(honeycombPayload, "[honeycomb.js] hcSendSpan()");
  event.add(honeycombPayload);
  event.send();
};

export const hcAddSpan = (data, stime) => {
  const duration = Date.now() - stime;
  const spanId = hcCreateId();
  const payload = {
    ...{ ...data },
    duration,
  }

  const trace = hcBuildTrace(hcRootTraceId, spanId, hcRootTraceId)
  hcSendSpan(payload, trace)
};

// export default {
//   rootTraceId,
//   honeycomb,
//   createId,
//   buildTrace,
//   sendSpan,
//   addSpan
// }