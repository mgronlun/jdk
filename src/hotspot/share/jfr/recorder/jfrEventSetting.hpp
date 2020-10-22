/*
 * Copyright (c) 2012, 2020, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 *
 */

#ifndef SHARE_JFR_RECORDER_JFREVENTSETTING_HPP
#define SHARE_JFR_RECORDER_JFREVENTSETTING_HPP

#include "jfr/utilities/jfrAllocation.hpp"
#include "jfrfiles/jfrEventControl.hpp"

//
// Native event settings as an associative array using the event id as key.
//
class JfrEventSetting : AllStatic {
 private:
  static JfrNativeSettings _jvm_event_settings;
  static jfrNativeEventSetting& setting(JfrEventId event_id);

 public:
  static void set_enabled(int64_t event_id, bool enabled);
  static bool is_enabled(JfrEventId event_id);
  static void set_stacktrace(int64_t event_id, bool enabled);
  static bool has_stacktrace(JfrEventId event_id);
  static bool set_threshold(int64_t event_id, int64_t threshold_ticks);
  static int64_t threshold(JfrEventId event_id);
  static bool set_cutoff(int64_t event_id, int64_t cutoff_ticks);
  static int64_t cutoff(JfrEventId event_id);
  static bool is_large(JfrEventId event_id);
  static void set_large(JfrEventId event_id);
  static bool set_ratelimit(int64_t event_id, int64_t ratelimit);
  static int64_t ratelimit(JfrEventId event_id);

  DEBUG_ONLY(static bool bounds_check_event(int64_t id);)
};

#endif // SHARE_JFR_RECORDER_JFREVENTSETTING_HPP
